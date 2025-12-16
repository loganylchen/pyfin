"""
POD5 file format parser
"""

from typing import List, Dict, Optional, Tuple, Any, Generator
from pathlib import Path
import logging

try:
    import pod5 as p5
except ImportError:
    raise ImportError("pod5 is required for POD5 support. " "Install it with: pip install pod5")


logger = logging.getLogger(__name__)


class Pod5Reader:
    """Reader for ONT POD5 files"""

    def __init__(self, file_path: str):
        """
        Initialize POD5 reader

        Args:
            file_path: Path to the pod5 file
        """
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"POD5 file not found: {file_path}")

        self._pod5_file = None
        self._read_ids = None
        self._read_count = None

    def __enter__(self):
        """Context manager entry"""
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.close()

    def open(self):
        """Open the pod5 file"""
        try:
            self._pod5_file = p5.Reader(self.file_path)
            self._read_ids = list(self._pod5_file.read_ids)
            self._read_count = len(self._read_ids)

            logger.info(f"Opened POD5 file with {self._read_count} reads: {self.file_path}")
        except Exception as e:
            logger.error(f"Failed to open POD5 file {self.file_path}: {e}")
            raise

    def close(self):
        """Close the pod5 file"""
        if self._pod5_file is not None:
            self._pod5_file.close()
            self._pod5_file = None
            logger.info(f"Closed POD5 file: {self.file_path}")

    @property
    def read_ids(self) -> List[str]:
        """
        Get list of all read IDs in the file

        Returns:
            List of read IDs
        """
        if self._read_ids is None:
            raise RuntimeError("POD5 file not opened. Call open() first.")
        return self._read_ids

    @property
    def read_count(self) -> int:
        """
        Get number of reads in the file

        Returns:
            Number of reads
        """
        if self._read_count is None:
            raise RuntimeError("POD5 file not opened. Call open() first.")
        return self._read_count

    def get_read(self, read_id: str) -> Optional[p5.Read]:
        """
        Get a Read object for a specific read ID

        Args:
            read_id: Read ID

        Returns:
            pod5.Read object or None if read not found
        """
        if self._pod5_file is None:
            raise RuntimeError("POD5 file not opened. Call open() first.")

        try:
            return self._pod5_file.reads(read_id)
        except KeyError:
            logger.warning(f"Read ID {read_id} not found in file")
            return None
        except Exception as e:
            logger.error(f"Failed to get read {read_id}: {e}")
            return None

    def get_all_read_data(self, read_id: str) -> Optional[Dict[str, Any]]:
        """
        Get all data for a specific read

        Args:
            read_id: Read ID

        Returns:
            Dictionary containing read data or None if read not found
        """
        read = self.get_read(read_id)
        if read is None:
            return None

        try:
            return {
                "read_id": read_id,
                "signal": read.signal.tolist(),
                "signal_shape": read.signal.shape,
                "signal_dtype": str(read.signal.dtype),
                "pore": read.pore,
                "calibration": {
                    "offset": float(read.calibration.offset),
                    "scale": float(read.calibration.scale),
                },
                "read_number": int(read.read_number),
                "start_sample": int(read.start_sample),
                "median_before": (
                    float(read.median_before) if read.median_before is not None else None
                ),
                "end_reason": int(read.end_reason) if read.end_reason is not None else None,
                "run_info": {
                    "acquisition_id": read.run_info.acquisition_id,
                    "acquisition_start_time": read.run_info.acquisition_start_time,
                    "adc_max": read.run_info.adc_max,
                    "adc_min": read.run_info.adc_min,
                    "context_tags": (
                        dict(read.run_info.context_tags) if read.run_info.context_tags else {}
                    ),
                    "device_id": read.run_info.device_id,
                    "digitisation": read.run_info.digitisation,
                    "flow_cell_id": read.run_info.flow_cell_id,
                    "flow_cell_product_code": read.run_info.flow_cell_product_code,
                    "protocol_name": read.run_info.protocol_name,
                    "protocol_run_id": read.run_info.protocol_run_id,
                    "protocol_start_time": read.run_info.protocol_start_time,
                    "sample_rate": read.run_info.sample_rate,
                    "sample_temperature": read.run_info.sample_temperature,
                    "sequencing_kit": read.run_info.sequencing_kit,
                    "track_codec": read.run_info.track_codec,
                    "version": read.run_info.version,
                    "well_name": (
                        read.run_info.well_name if hasattr(read.run_info, "well_name") else None
                    ),
                },
                "tracked_scaling": (
                    {
                        "scale": (
                            float(read.tracked_scaling.scale) if read.tracked_scaling else None
                        ),
                        "shift": (
                            float(read.tracked_scaling.shift) if read.tracked_scaling else None
                        ),
                        "digitisation": (
                            read.tracked_scaling.digitisation if read.tracked_scaling else None
                        ),
                        "variance": (
                            float(read.tracked_scaling.variance) if read.tracked_scaling else None
                        ),
                    }
                    if read.tracked_scaling
                    else None
                ),
                "num_minknow_events": (
                    read.num_minknow_events if hasattr(read, "num_minknow_events") else None
                ),
            }

        except Exception as e:
            logger.error(f"Failed to extract data from read {read_id}: {e}")
            return None

    def get_signal(self, read_id: str) -> Optional[Tuple[List[int], Dict[str, Any]]]:
        """
        Get raw signal for a specific read

        Args:
            read_id: Read ID

        Returns:
            Tuple of (signal_array, metadata_dict) or None if read not found
        """
        read = self.get_read(read_id)
        if read is None:
            return None

        metadata = {
            "read_id": read_id,
            "duration": len(read.signal),
            "sample_rate": read.run_info.sample_rate,
            "channel_id": read.pore,
            "start_sample": read.start_sample,
            "read_number": read.read_number,
            "calibration": {
                "offset": float(read.calibration.offset),
                "scale": float(read.calibration.scale),
            },
            "digitisation": read.run_info.digitisation,
            "median_before": float(read.median_before) if read.median_before is not None else None,
        }

        return read.signal.tolist(), metadata

    def get_calibrated_signal(self, read_id: str) -> Optional[Tuple[List[float], Dict[str, Any]]]:
        """
        Get calibrated signal (converted to picoamperes) for a specific read

        Args:
            read_id: Read ID

        Returns:
            Tuple of (picoamp_array, metadata_dict) or None if read not found
        """
        read = self.get_read(read_id)
        if read is None:
            return None

        try:
            # Apply calibration to convert to picoamperes
            calibrated_signal = read.signal.astype(float)
            calibrated_signal = (
                calibrated_signal + read.calibration.offset
            ) * read.calibration.scale

            metadata = {
                "read_id": read_id,
                "duration": len(read.signal),
                "sample_rate": read.run_info.sample_rate,
                "channel_id": read.pore,
                "start_sample": read.start_sample,
                "read_number": read.read_number,
                "calibration": {
                    "offset": float(read.calibration.offset),
                    "scale": float(read.calibration.scale),
                },
                "median_before": (
                    float(read.median_before) if read.median_before is not None else None
                ),
                "unit": "picoamperes",
            }

            return calibrated_signal.tolist(), metadata

        except Exception as e:
            logger.error(f"Failed to calibrate signal for read {read_id}: {e}")
            return None

    def get_batch_reads(self, read_ids: List[str]) -> List[Optional[p5.Read]]:
        """
        Get multiple reads in batch (more efficient than individual reads)

        Args:
            read_ids: List of read IDs to retrieve

        Returns:
            List of pod5.Read objects or None for each read ID
        """
        if self._pod5_file is None:
            raise RuntimeError("POD5 file not opened. Call open() first.")

        try:
            return self._pod5_file.reads(read_ids)
        except Exception as e:
            logger.error(f"Failed to get batch reads: {e}")
            return [None] * len(read_ids)

    def iterate_reads(self, batch_size: int = 1000) -> Generator[Dict[str, Any], None, None]:
        """
        Generator to iterate through all reads in the file

        Args:
            batch_size: Number of reads to fetch per batch (default: 1000)

        Yields:
            Read data dictionary for each read
        """
        if self._read_ids is None:
            raise RuntimeError("POD5 file not opened. Call open() first.")

        for i in range(0, len(self._read_ids), batch_size):
            batch = self._read_ids[i : i + batch_size]
            reads = self.get_batch_reads(batch)

            for read in reads:
                if read is not None:
                    read_data = self.get_all_read_data(read.read_id)
                    if read_data:
                        yield read_data

    def get_file_info(self) -> Dict[str, Any]:
        """
        Get file information

        Returns:
            Dictionary with file information
        """
        if self._pod5_file is None:
            raise RuntimeError("POD5 file not opened. Call open() first.")

        try:
            run_info = self._pod5_file.run_info
            return {
                "filename": str(self.file_path),
                "num_reads": self._read_count,
                "num_chunks": self._pod5_file.num_chunks,
                "compression": str(self._pod5_file.compression),
                "version": self._pod5_file.version,
                "run_info": {
                    "flow_cell_id": run_info.flow_cell_id,
                    "sample_id": run_info.sample_id if hasattr(run_info, "sample_id") else None,
                    "sequencing_kit": run_info.sequencing_kit,
                    "protocol_name": run_info.protocol_name,
                    "device_id": run_info.device_id,
                    "sample_rate": run_info.sample_rate,
                    "acquisition_id": run_info.acquisition_id,
                    "acquisition_start_time": run_info.acquisition_start_time,
                },
            }
        except Exception as e:
            logger.error(f"Failed to get file info: {e}")
            return {}

    @staticmethod
    def convert_fast5_to_pod5(
        fast5_paths: List[str], output_path: str, recursive: bool = False
    ) -> bool:
        """
        Convert FAST5 files to POD5 (requires pod5 command line tools)

        Args:
            fast5_paths: List of paths to FAST5 files or directories
            output_path: Output POD5 file path
            recursive: Search directories recursively

        Returns:
            True if conversion successful, False otherwise
        """
        try:
            import subprocess

            cmd = ["pod5", "convert", "fast5"]
            for path in fast5_paths:
                cmd.extend(["--input", path])
            cmd.extend(["--output", output_path])
            if recursive:
                cmd.append("--recursive")

            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                logger.info(f"Successfully converted FAST5 to POD5: {output_path}")
                return True
            else:
                logger.error(f"Failed to convert: {result.stderr}")
                return False
        except FileNotFoundError:
            logger.error("pod5 command not found. Install pod5 tools for format conversion.")
            return False
        except Exception as e:
            logger.error(f"Error during conversion: {e}")
            return False
