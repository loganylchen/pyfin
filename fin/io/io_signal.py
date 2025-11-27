"""
Signal file parsers for fast5, pod5, and blow5/slow5 formats.

Handles reading raw nanopore ionic current signals from various file formats,
extracting both signal data and associated metadata.
"""

import h5py
import pod5
import numpy as np
from typing import Dict, List, Optional, Iterator, Union, Any, Tuple
from pathlib import Path
import logging
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class ReadSignal:
    """Container for a read's raw signal data."""
    read_id: str
    signal: np.ndarray
    sampling_rate: float
    start_time: float
    read_number: int
    channel_number: int
    mux: int
    median_before: float
    pore_type: Optional[str] = None
    sequence_length: Optional[int] = None
    context_tags: Optional[Dict[str, str]] = None
    tracking_id: Optional[Dict[str, str]] = None


def is_fast5_file(filepath: str) -> bool:
    """Check if file is a fast5 file."""
    return str(filepath).lower().endswith('.fast5')


def is_multi_fast5_file(filepath: str) -> bool:
    """Check if file is a multi-read fast5 file."""
    if not is_fast5_file(filepath):
        return False

    try:
        with h5py.File(filepath, 'r') as f:
            # Multi-fast5 has 'read_x' groups at root level
            for key in f.keys():
                if key.startswith('read_'):
                    return True
            return False
    except OSError:
        return False


def is_single_fast5_file(filepath: str) -> bool:
    """Check if file is a single-read fast5 file."""
    if not is_fast5_file(filepath):
        return False

    try:
        with h5py.File(filepath, 'r') as f:
            # Single fast5 has Raw/Reads/Read_X pattern
            if 'Raw' in f and 'Reads' in f['Raw']:
                return True
            return False
    except OSError:
        return False


def is_pod5_file(filepath: str) -> bool:
    """Check if file is a pod5 file."""
    return str(filepath).lower().endswith('.pod5')


def is_slow5_file(filepath: str) -> bool:
    """Check if file has slow5/blow5 extension."""
    filepath = str(filepath).lower()
    return filepath.endswith('.slow5') or filepath.endswith('.blow5')


class Fast5Parser:
    """Parser for Oxford Nanopore fast5 files (HDF5 format)."""

    def __init__(self, filepath: str):
        """
        Initialize fast5 parser.

        Args:
            filepath: Path to fast5 file (single or multi-read)
        """
        self.filepath = filepath
        self.is_multi = is_multi_fast5_file(filepath)
        self._file_handle = None

        logger.info(f"Created Fast5Parser for {filepath} (multi={self.is_multi})")

    def __enter__(self):
        """Context manager entry."""
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()

    def open(self):
        """Open the fast5 file."""
        if self._file_handle is not None:
            return

        try:
            self._file_handle = h5py.File(self.filepath, 'r')
            logger.debug(f"Opened fast5 file: {self.filepath}")
        except OSError as e:
            logger.error(f"Failed to open fast5 file {self.filepath}: {e}")
            raise

    def close(self):
        """Close the fast5 file."""
        if self._file_handle is not None:
            self._file_handle.close()
            self._file_handle = None
            logger.debug(f"Closed fast5 file: {self.filepath}")

    def get_read_ids(self) -> List[str]:
        """Get list of all read IDs in the file."""
        if self._file_handle is None:
            raise ValueError("Fast5 file not open")

        read_ids = []

        if self.is_multi:
            # Multi-read fast5: read groups at root level
            for group_name in self._file_handle.keys():
                if group_name.startswith('read_'):
                    # Extract read_id from group attributes
                    group = self._file_handle[group_name]
                    if 'read_id' in group.attrs:
                        read_ids.append(group.attrs['read_id'])
                    else:
                        logger.warning(f"No read_id in group {group_name}")
        else:
            # Single-read fast5: one read in Raw/Reads/
            if 'Raw' in self._file_handle and 'Reads' in self._file_handle['Raw']:
                reads_group = self._file_handle['Raw']['Reads']
                for read_name in reads_group.keys():
                    if 'read_id' in reads_group[read_name].attrs:
                        read_ids.append(reads_group[read_name].attrs['read_id'])

        logger.info(f"Found {len(read_ids)} reads in {self.filepath}")
        return read_ids

    def _extract_signal_single_read(
        self,
        read_group: h5py.Group
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """
        Extract signal from single read group.

        Args:
            read_group: HDF5 group containing read data

        Returns:
            Tuple of (signal_array, metadata_dict)
        """
        # Get signal dataset
        signal_ds = read_group['Signal']
        signal = np.array(signal_ds, dtype=np.int16)

        # Get scaling parameters
        if 'range' in read_group['Signal'].attrs:
            # Newer format: range and digitisation
            range_val = read_group['Signal'].attrs['range']
            digitisation = read_group['Signal'].attrs['digitisation']
            offset = read_group['Signal'].attrs.get('offset', 0)

            # Convert from raw ADC to pA
            signal_pa = signal * (range_val / digitisation) + offset
        else:
            # Older format: try to find scaling in attributes
            signal_pa = signal.astype(np.float32)
            logger.warning("Could not find signal scaling parameters")

        # Extract metadata
        metadata = {}

        # Read attributes
        for attr_name, attr_value in read_group.attrs.items():
            if isinstance(attr_value, (str, int, float)):
                metadata[attr_name] = attr_value

        # Read channel_id information
        if 'channel_id' in read_group:
            channel_id = read_group['channel_id']
            for attr_name, attr_value in channel_id.attrs.items():
                if hasattr(attr_value, 'decode'):
                    metadata[attr_name] = attr_value.decode('utf-8')
                else:
                    metadata[attr_name] = attr_value

        return signal_pa, metadata

    def get_read_signal(self, read_id: str) -> Optional[ReadSignal]:
        """
        Extract signal for a specific read ID.

        Args:
            read_id: Read identifier

        Returns:
            ReadSignal object or None if not found
        """
        if self._file_handle is None:
            raise ValueError("Fast5 file not open")

        found = False
        signal = None
        metadata = {}

        if self.is_multi:
            # Multi-read: search for read group with matching ID
            for group_name in self._file_handle.keys():
                if not group_name.startswith('read_'):
                    continue

                group = self._file_handle[group_name]
                if 'read_id' in group.attrs and group.attrs['read_id'] == read_id:
                    found = True

                    # Extract signal from Raw subgroup
                    if 'Raw' in group and 'Signal' in group['Raw']:
                        read_group = group['Raw']
                        signal, md = self._extract_signal_single_read(read_group)
                        metadata.update(md)
                    break
        else:
            # Single-read: check if this is the right read
            if 'Raw' in self._file_handle and 'Reads' in self._file_handle['Raw']:
                reads_group = self._file_handle['Raw']['Reads']
                for read_name in reads_group.keys():
                    read_group = reads_group[read_name]
                    if 'read_id' in read_group.attrs and read_group.attrs['read_id'] == read_id:
                        found = True
                        signal, md = self._extract_signal_single_read(read_group)
                        metadata.update(md)
                        break

        if not found:
            logger.error(f"Read ID {read_id} not found in {self.filepath}")
            return None

        if signal is None:
            logger.error(f"Failed to extract signal for {read_id}")
            return None

        # Create ReadSignal object
        read_signal = ReadSignal(
            read_id=read_id,
            signal=signal,
            sampling_rate=metadata.get('sampling_rate', 4000.0),
            start_time=metadata.get('start_time', 0.0),
            read_number=metadata.get('read_number', 0),
            channel_number=metadata.get('channel_number', 0),
            mux=metadata.get('mux', 0),
            median_before=metadata.get('median_before', 0.0),
            pore_type=metadata.get('pore_type', None),
            sequence_length=metadata.get('sequence_length', None),
            context_tags=self._extract_context_tags(),
            tracking_id=self._extract_tracking_id()
        )

        return read_signal

    def _extract_context_tags(self) -> Dict[str, str]:
        """Extract context tags from file."""
        if self._file_handle is None:
            return {}

        context_tags = {}
        try:
            if 'UniqueGlobalKey/context_tags' in self._file_handle:
                tags_group = self._file_handle['UniqueGlobalKey/context_tags']
                for key, value in tags_group.attrs.items():
                    if hasattr(value, 'decode'):
                        context_tags[key] = value.decode('utf-8')
                    else:
                        context_tags[key] = str(value)
        except Exception as e:
            logger.debug(f"Could not extract context tags: {e}")

        return context_tags

    def _extract_tracking_id(self) -> Dict[str, str]:
        """Extract tracking ID from file."""
        if self._file_handle is None:
            return {}

        tracking_id = {}
        try:
            if 'UniqueGlobalKey/tracking_id' in self._file_handle:
                track_group = self._file_handle['UniqueGlobalKey/tracking_id']
                for key, value in track_group.attrs.items():
                    if hasattr(value, 'decode'):
                        tracking_id[key] = value.decode('utf-8')
                    else:
                        tracking_id[key] = str(value)
        except Exception as e:
            logger.debug(f"Could not extract tracking ID: {e}")

        return tracking_id

    def iterate_reads(self) -> Iterator[ReadSignal]:
        """
        Iterate over all reads in the file.

        Yields:
            ReadSignal objects
        """
        read_ids = self.get_read_ids()

        for read_id in read_ids:
            read_signal = self.get_read_signal(read_id)
            if read_signal is not None:
                yield read_signal


class Pod5Parser:
    """Parser for ONT Pod5 files (next-gen format replacing fast5)."""

    def __init__(self, filepath: str):
        """
        Initialize Pod5 parser.

        Args:
            filepath: Path to Pod5 file
        """
        self.filepath = filepath
        self.reader = None

    def __enter__(self):
        """Context manager entry."""
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()

    def open(self):
        """Open the Pod5 file."""
        if self.reader is not None:
            return

        try:
            self.reader = pod5.Reader(self.filepath)
            logger.info(f"Opened Pod5 file: {self.filepath}")
        except Exception as e:
            logger.error(f"Failed to open Pod5 file {self.filepath}: {e}")
            raise

    def close(self):
        """Close the Pod5 file."""
        if self.reader is not None:
            self.reader = None
            logger.debug(f"Closed Pod5 file: {self.filepath}")

    def get_read_ids(self) -> List[str]:
        """Get list of all read IDs."""
        if self.reader is None:
            raise ValueError("Pod5 file not open")

        read_ids = [str(read.read_id) for read in self.reader.reads()]
        logger.info(f"Found {len(read_ids)} reads in {self.filepath}")
        return read_ids

    def get_read_signal(self, read_id: str) -> Optional[ReadSignal]:
        """
        Extract signal for specific read ID.

        Args:
            read_id: Read identifier

        Returns:
            ReadSignal object
        """
        if self.reader is None:
            raise ValueError("Pod5 file not open")

        try:
            read = self.reader.get_read(read_id)

            # Convert signal to pA
            signal_pA = read.signal_as_pA()

            # Create ReadSignal
            read_signal = ReadSignal(
                read_id=str(read.read_id),
                signal=signal_pA,
                sampling_rate=read.run_info.sample_rate,
                start_time=read.start_time,
                read_number=read.read_number,
                channel_number=read.pore.channel,
                mux=read.pore.well,
                median_before=read.median_before,
                pore_type=read.pore.pore_type,
                sequence_length=read.num_samples,
                context_tags={},
                tracking_id={}
            )

            return read_signal

        except Exception as e:
            logger.error(f"Failed to extract read {read_id}: {e}")
            return None

    def iterate_reads(self) -> Iterator[ReadSignal]:
        """Iterate over all reads."""
        if self.reader is None:
            raise ValueError("Pod5 file not open")

        for read_record in self.reader.reads():
            read_signal = self.get_read_signal(str(read_record.read_id))
            if read_signal is not None:
                yield read_signal


def get_signal_parser(filepath: str):
    """
    Factory function to get appropriate parser for file type.

    Args:
        filepath: Path to signal file

    Returns:
        Parser object (Fast5Parser or Pod5Parser)
    """
    if is_fast5_file(filepath):
        return Fast5Parser(filepath)
    elif is_pod5_file(filepath):
        return Pod5Parser(filepath)
    else:
        raise ValueError(f"Unsupported file format: {filepath}")


def get_read_signal_from_file(filepath: str, read_id: str) -> Optional[ReadSignal]:
    """
    Convenience function to extract signal for a specific read.

    Args:
        filepath: Path to signal file
        read_id: Read identifier

    Returns:
        ReadSignal object or None
    """
    parser_class = get_signal_parser(filepath)

    with parser_class as parser:
        return parser.get_read_signal(read_id)


def iterate_reads_from_file(filepath: str) -> Iterator[ReadSignal]:
    """
    Convenience function to iterate over all reads in a file.

    Args:
        filepath: Path to signal file

    Yields:
        ReadSignal objects
    """
    parser_class = get_signal_parser(filepath)

    with parser_class as parser:
        yield from parser.iterate_reads()
