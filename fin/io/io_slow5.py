"""
Slow5 file format parser using pyslow5
"""

from typing import List, Dict, Optional, Tuple, Any, Generator
from pathlib import Path
import logging

try:
    import pyslow5
except ImportError:
    raise ImportError(
        "pyslow5 is required for Slow5 support. "
        "Install it with: pip install pyslow5"
    )


logger = logging.getLogger(__name__)


class Slow5Reader:
    """Reader for SLOW5/BLOW5 files"""

    def __init__(self, file_path: str):
        """
        Initialize Slow5 reader

        Args:
            file_path: Path to the slow5/blow5 file
        """
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"Slow5 file not found: {file_path}")

        self._s5_file = None
        self._header = None
        self._read_count = None

    def __enter__(self):
        """Context manager entry"""
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.close()

    def open(self):
        """Open the slow5 file"""
        try:
            self._s5_file = pyslow5.Open(str(self.file_path), 'r')
            self._header = self._s5_file.get_all_headers()

            # Count reads
            read_ids = self._s5_file.get_read_ids()
            self._read_count = len(read_ids) if read_ids else 0

            logger.info(f"Opened slow5 file with {self._read_count} reads: {self.file_path}")
        except Exception as e:
            logger.error(f"Failed to open slow5 file {self.file_path}: {e}")
            raise

    def close(self):
        """Close the slow5 file"""
        if self._s5_file is not None:
            self._s5_file.close()
            self._s5_file = None
            logger.info(f"Closed slow5 file: {self.file_path}")

    @property
    def read_ids(self) -> List[str]:
        """
        Get list of all read IDs in the file

        Returns:
            List of read IDs
        """
        if self._s5_file is None:
            raise RuntimeError("Slow5 file not opened. Call open() first.")

        try:
            return self._s5_file.get_read_ids()
        except Exception as e:
            logger.error(f"Failed to get read IDs: {e}")
            return []

    @property
    def header(self) -> Dict[str, Any]:
        """
        Get file header information

        Returns:
            Dictionary containing header information
        """
        if self._header is None:
            raise RuntimeError("Slow5 file not opened. Call open() first.")
        return dict(self._header)

    def get_read(self, read_id: str) -> Optional[Dict[str, Any]]:
        """
        Get all data for a specific read

        Args:
            read_id: Read ID

        Returns:
            Dictionary containing read data or None if read not found
        """
        if self._s5_file is None:
            raise RuntimeError("Slow5 file not opened. Call open() first.")

        try:
            # Get read data (returns numpy arrays)
            read_data = self._s5_file.get_read(read_id)

            if read_data is None:
                logger.warning(f"Read ID {read_id} not found in file")
                return None

            return {
                'read_id': read_id,
                'signal': read_data['signal'].tolist(),
                'raw_signal': read_data['raw_signal'].tolist() if 'raw_signal' in read_data else None,
                'picoamp_values': read_data['picoamp_values'].tolist() if 'picoamp_values' in read_data else None,
                'digitisation': read_data['digitisation'],
                'offset': read_data['offset'],
                'range': read_data['range'],
                'channel': read_data.get('channel_number') or read_data.get('channel'),
                'sampling_rate': read_data['sampling_rate'],
                'start_time': read_data['start_time'],
                'read_number': read_data['read_number'],
                'start_mux': read_data['start_mux'],
                'median_before': read_data.get('median_before', None),
                'median_after': read_data.get('median_after', None),
            }

        except Exception as e:
            logger.warning(f"Failed to get read {read_id}: {e}")
            return None

    def get_signal(self, read_id: str) -> Optional[Tuple[List[int], Dict[str, Any]]]:
        """
        Get raw signal for a specific read

        Args:
            read_id: Read ID

        Returns:
            Tuple of (signal_array, metadata_dict) or None if read not found
        """
        read_data = self.get_read(read_id)
        if read_data is None:
            return None

        metadata = {
            'read_id': read_id,
            'duration': len(read_data['signal']),
            'sample_rate': read_data['sampling_rate'],
            'channel_id': read_data['channel'],
            'start_time': read_data['start_time'],
            'digitisation': read_data['digitisation'],
            'offset': read_data['offset'],
            'range': read_data['range'],
        }

        return read_data['signal'], metadata

    def get_raw_signal(self, read_id: str) -> Optional[Tuple[List[int], Dict[str, Any]]]:
        """
        Get raw signal values (before scaling) for a specific read

        Args:
            read_id: Read ID

        Returns:
            Tuple of (raw_signal_array, metadata_dict) or None if read not found
        """
        read_data = self.get_read(read_id)
        if read_data is None or read_data['raw_signal'] is None:
            return None

        metadata = {
            'read_id': read_id,
            'duration': len(read_data['raw_signal']),
            'sample_rate': read_data['sampling_rate'],
            'channel_id': read_data['channel'],
            'start_time': read_data['start_time'],
        }

        return read_data['raw_signal'], metadata

    def get_picoamp_signal(self, read_id: str) -> Optional[Tuple[List[float], Dict[str, Any]]]:
        """
        Get picoampere-scaled signal for a specific read

        Args:
            read_id: Read ID

        Returns:
            Tuple of (picoamp_array, metadata_dict) or None if read not found
        """
        read_data = self.get_read(read_id)
        if read_data is None or read_data['picoamp_values'] is None:
            return None

        metadata = {
            'read_id': read_id,
            'duration': len(read_data['picoamp_values']),
            'sample_rate': read_data['sampling_rate'],
            'channel_id': read_data['channel'],
            'start_time': read_data['start_time'],
            'digitisation': read_data['digitisation'],
            'offset': read_data['offset'],
            'range': read_data['range'],
        }

        return read_data['picoamp_values'], metadata

    def iterate_reads(self) -> Generator[Dict[str, Any], None, None]:
        """
        Generator to iterate through all reads in the file

        Yields:
            Read data dictionary for each read
        """
        if self._s5_file is None:
            raise RuntimeError("Slow5 file not opened. Call open() first.")

        read_ids = self.read_ids
        if not read_ids:
            logger.warning("No reads found in file")
            return

        for read_id in read_ids:
            read_data = self.get_read(read_id)
            if read_data:
                yield read_data

    def get_file_info(self) -> Dict[str, Any]:
        """
        Get file information

        Returns:
            Dictionary with file information
        """
        if self._header is None:
            raise RuntimeError("Slow5 file not opened. Call open() first.")

        try:
            return {
                'filename': str(self.file_path),
                'num_reads': self._read_count,
                'file_version': self._header.get('file_version', 'unknown'),
                'format': 'SLOW5' if self.file_path.suffix.lower() == '.slow5' else 'BLOW5',
                'header': dict(self._header)
            }
        except Exception as e:
            logger.error(f"Failed to get file info: {e}")
            return {}

    @staticmethod
    def convert_to_slow5(blow5_path: str, slow5_path: str) -> bool:
        """
        Convert BLOW5 file to SLOW5 (requires slow5tools)

        Args:
            blow5_path: Path to input BLOW5 file
            slow5_path: Path to output SLOW5 file

        Returns:
            True if conversion successful, False otherwise
        """
        try:
            # This uses slow5tools command line if available
            import subprocess
            result = subprocess.run(
                ['slow5tools', 'view', blow5_path, '-o', slow5_path],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                logger.info(f"Successfully converted {blow5_path} to {slow5_path}")
                return True
            else:
                logger.error(f"Failed to convert: {result.stderr}")
                return False
        except FileNotFoundError:
            logger.error("slow5tools not found. Install it for format conversion.")
            return False
        except Exception as e:
            logger.error(f"Error during conversion: {e}")
            return False
