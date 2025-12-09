"""
Fast5 file format parser using ont-fast5-api
"""

from typing import List, Dict, Optional, Tuple, Any
from pathlib import Path
import logging

try:
    from ont_fast5_api.fast5_interface import get_fast5_file
    from ont_fast5_api.analysis_tools import basecall_1d
except ImportError:
    raise ImportError(
        "ont-fast5-api is required for Fast5 support. "
        "Install it with: pip install ont-fast5-api"
    )


logger = logging.getLogger(__name__)


class Fast5Reader:
    """Reader for Oxford Nanopore FAST5 files"""

    def __init__(self, file_path: str):
        """
        Initialize Fast5 reader

        Args:
            file_path: Path to the fast5 file
        """
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"Fast5 file not found: {file_path}")

        self._fast5_file = None
        self._read_ids = None

    def __enter__(self):
        """Context manager entry"""
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.close()

    def open(self):
        """Open the fast5 file"""
        try:
            self._fast5_file = get_fast5_file(self.file_path, mode='r')
            self._read_ids = list(self._fast5_file.get_read_ids())
            logger.info(f"Opened fast5 file with {len(self._read_ids)} reads: {self.file_path}")
        except Exception as e:
            logger.error(f"Failed to open fast5 file {self.file_path}: {e}")
            raise

    def close(self):
        """Close the fast5 file"""
        if self._fast5_file is not None:
            self._fast5_file.close()
            self._fast5_file = None
            logger.info(f"Closed fast5 file: {self.file_path}")

    @property
    def read_ids(self) -> List[str]:
        """
        Get list of all read IDs in the file

        Returns:
            List of read IDs
        """
        if self._read_ids is None:
            raise RuntimeError("Fast5 file not opened. Call open() first.")
        return self._read_ids

    def get_raw_signal(self, read_id: str) -> Optional[Tuple[List[int], Dict[str, Any]]]:
        """
        Get raw signal for a specific read

        Args:
            read_id: Read ID

        Returns:
            Tuple of (signal_array, metadata_dict) or None if read not found
        """
        if self._fast5_file is None:
            raise RuntimeError("Fast5 file not opened. Call open() first.")

        try:
            read = self._fast5_file.get_read(read_id)

            # Get raw data
            raw_data = read.get_raw_data()

            # Get channel info
            channel_info = read.get_channel_info()

            metadata = {
                'read_id': read_id,
                'duration': len(raw_data),
                'sample_rate': channel_info.get('sampling_rate'),
                'channel_id': channel_info.get('channel_number'),
                'start_time': channel_info.get('start_time'),
            }

            return list(raw_data), metadata

        except Exception as e:
            logger.warning(f"Failed to get raw signal for read {read_id}: {e}")
            return None

    def get_basecall_data(self, read_id: str, analysis_group: str = "Basecall_1D_000") -> Optional[Dict[str, Any]]:
        """
        Get basecall data for a specific read

        Args:
            read_id: Read ID
            analysis_group: Analysis group name (default: "Basecall_1D_000")

        Returns:
            Dictionary containing basecall data or None if not found
        """
        if self._fast5_file is None:
            raise RuntimeError("Fast5 file not opened. Call open() first.")

        try:
            read = self._fast5_file.get_read(read_id)

            # Get basecall data
            basecall_data = read.get_analysis_dataset(analysis_group, 'BaseCalled_template')
            fastq = read.get_analysis_dataset(analysis_group, 'BaseCalled_template/Fastq')

            # Get move table if available
            try:
                move_table = read.get_analysis_dataset(analysis_group, 'BaseCalled_template/Move')
            except:
                move_table = None

            # Get trace if available
            try:
                trace = read.get_analysis_dataset(analysis_group, 'BaseCalled_template/Trace')
            except:
                trace = None

            return {
                'fastq': fastq,
                'basecall_data': basecall_data,
                'move_table': move_table,
                'trace': trace,
                'analysis_group': analysis_group
            }

        except Exception as e:
            logger.warning(f"Failed to get basecall data for read {read_id}: {e}")
            return None

    def get_all_read_data(self, read_id: str) -> Optional[Dict[str, Any]]:
        """
        Get all available data for a specific read

        Args:
            read_id: Read ID

        Returns:
            Dictionary containing all read data or None if read not found
        """
        result = {
            'read_id': read_id,
            'raw_signal': None,
            'basecall': None
        }

        raw_data = self.get_raw_signal(read_id)
        if raw_data:
            result['raw_signal'] = raw_data[0]
            result['raw_metadata'] = raw_data[1]

        basecall_data = self.get_basecall_data(read_id)
        if basecall_data:
            result['basecall'] = basecall_data

        return result if result['raw_signal'] or result['basecall'] else None

    def iterate_reads(self):
        """
        Generator to iterate through all reads in the file

        Yields:
            Read data dictionary for each read
        """
        if self._read_ids is None:
            raise RuntimeError("Fast5 file not opened. Call open() first.")

        for read_id in self._read_ids:
            read_data = self.get_all_read_data(read_id)
            if read_data:
                yield read_data

    @staticmethod
    def get_version_info(file_path: str) -> Dict[str, str]:
        """
        Get version information from a fast5 file

        Args:
            file_path: Path to the fast5 file

        Returns:
            Dictionary with version info
        """
        try:
            with get_fast5_file(file_path, mode='r') as f5:
                version_info = f5.status.get('ont-fast5-api', {})
                return {
                    'file_version': version_info.get('file_format_version', 'unknown'),
                    'api_version': version_info.get('api_version', 'unknown'),
                    'tracking_id': dict(f5.status.get('tracking_id', {}))
                }
        except Exception as e:
            logger.error(f"Failed to get version info for {file_path}: {e}")
            return {}
