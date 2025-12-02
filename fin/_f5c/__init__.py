"""
f5c Eventalign Module

This module provides Python bindings for f5c's eventalign functionality,
which aligns nanopore signal events to reference k-mers.
"""

# from typing import List, Dict, Iterator, Optional, Tuple
# from pathlib import Path
# import warnings

# from ..utils.log_config import get_package_logger

# logger = get_package_logger(__name__)

# try:
#     from . import _f5c
#     _F5C_AVAILABLE = True
# except ImportError:
#     _F5C_AVAILABLE = False
#     logger.warning("f5c module not available. Install with: pip install -e .")


# class EventAligner:
#     """
#     Python wrapper for f5c eventalign functionality.

#     This class provides a clean interface to f5c's eventalign function,
#     which aligns nanopore signal events to reference k-mers.

#     Example:
#         >>> aligner = EventAligner(
#         ...     bam_path="reads.bam",
#         ...     fasta_path="reference.fa",
#         ...     slow5_path="signals.blow5"
#         ... )
#         >>> for alignment in aligner.align_reads():
#         ...     print(f"Read: {alignment['read_id']}")
#         ...     for event in alignment['events']:
#         ...         print(f"  Kmer: {event['kmer']}, Event mean: {event['event_mean']}")
#     """

#     def __init__(self, bam_path: str, fasta_path: str, slow5_path: Optional[str] = None):
#         """
#         Initialize EventAligner.

#         Args:
#             bam_path: Path to BAM file with read alignments
#             fasta_path: Path to FASTA reference file
#             slow5_path: Optional path to SLOW5/BLOW5 signal file

#         Raises:
#             ImportError: If f5c module is not compiled
#             FileNotFoundError: If input files don't exist
#         """
#         if not _F5C_AVAILABLE:
#             raise ImportError(
#                 "f5c module not available. Please compile with: pip install -e ."
#             )

#         # Validate files
#         for path in [bam_path, fasta_path]:
#             if not Path(path).exists():
#                 raise FileNotFoundError(f"File not found: {path}")

#         if slow5_path and not Path(slow5_path).exists():
#             raise FileNotFoundError(f"File not found: {slow5_path}")

#         self.bam_path = str(Path(bam_path).resolve())
#         self.fasta_path = str(Path(fasta_path).resolve())
#         self.slow5_path = str(Path(slow5_path).resolve()) if slow5_path else None

#         # Initialize f5c context
#         try:
#             self._context = _f5c.init_eventalign(
#                 self.bam_path,
#                 self.fasta_path,
#                 self.slow5_path
#             )
#             logger.info(f"Initialized EventAligner with BAM: {bam_path}")
#         except Exception as e:
#             logger.error(f"Failed to initialize EventAligner: {e}")
#             raise

#     def align_batch(self) -> List[Dict]:
#         """
#         Align events for the next batch of reads.

#         Returns:
#             List of alignment dictionaries, each containing:
#             - event_idx: Index of the event
#             - kmer_idx: Index of the k-mer
#             - kmer: K-mer sequence
#             - event_mean: Mean of the event signal
#             - event_stdv: Standard deviation of the event signal
#             - model_mean: Expected mean from the model
#             - model_stdv: Expected stdv from the model
#             - posterior_probability: Posterior probability of the alignment
#             - start_idx: Start index in the signal
#             - end_idx: End index in the signal

#         Example:
#             >>> alignments = aligner.align_batch()
#             >>> for aln in alignments:
#             ...     print(f"{aln['kmer']}: {aln['event_mean']:.2f}")
#         """
#         if not self._context:
#             raise RuntimeError("EventAligner not initialized")

#         try:
#             results = _f5c.eventalign_read(self._context, "")
#             return results
#         except Exception as e:
#             logger.error(f"Error during event alignment: {e}")
#             raise

#     def align_reads(self, max_reads: Optional[int] = None) -> Iterator[Tuple[str, List[Dict]]]:
#         """
#         Generator that yields event alignments for all reads.

#         Args:
#             max_reads: Optional maximum number of reads to process

#         Yields:
#             Tuple of (read_id, alignments) for each read

#         Example:
#             >>> for read_id, alignments in aligner.align_reads():
#             ...     print(f"Read {read_id}: {len(alignments)} alignments")
#             ...     for aln in alignments[:5]:
#             ...         print(f"  {aln['kmer']}: {aln['event_mean']:.2f}")
#         """
#         read_count = 0

#         while True:
#             if max_reads is not None and read_count >= max_reads:
#                 break

#             try:
#                 alignments = self.align_batch()

#                 if not alignments:
#                     break

#                 # Group by read (simplified - assumes all alignments in batch belong to same read)
#                 # In real implementation, you might want to track read IDs from the BAM
#                 read_id = f"read_{read_count}"
#                 yield read_id, alignments

#                 read_count += 1

#                 if read_count % 100 == 0:
#                     logger.info(f"Processed {read_count} reads")

#             except Exception as e:
#                 logger.error(f"Error processing read {read_count}: {e}")
#                 break

#     def get_summary_stats(self) -> Dict:
#         """
#         Get summary statistics for the alignments processed so far.

#         Returns:
#             Dictionary with:
#             - total_reads: Number of reads processed
#             - total_events: Total number of events aligned
#             - avg_events_per_read: Average events per read

#         Note:
#             This is a placeholder. In a full implementation, you would
#             track these statistics during alignment.
#         """
#         return {
#             "total_reads": 0,
#             "total_events": 0,
#             "avg_events_per_read": 0.0,
#         }

#     def __enter__(self):
#         """Context manager entry"""
#         return self

#     def __exit__(self, exc_type, exc_val, exc_tb):
#         """Context manager exit"""
#         self.close()

#     def close(self):
#         """Release resources"""
#         if hasattr(self, '_context') and self._context:
#             try:
#                 _f5c.free_eventalign(self._context)
#                 self._context = None
#                 logger.info("EventAligner resources released")
#             except Exception as e:
#                 logger.warning(f"Error releasing resources: {e}")

#     def __del__(self):
#         """Destructor"""
#         self.close()


# def eventalign(bam_path: str,
#                fasta_path: str,
#                slow5_path: Optional[str] = None,
#                max_reads: Optional[int] = None) -> Iterator[Tuple[str, List[Dict]]]:
#     """
#     Convenience function for event alignment.

#     This is a generator function that creates an EventAligner and yields
#     alignments for all reads.

#     Args:
#         bam_path: Path to BAM file
#         fasta_path: Path to FASTA reference
#         slow5_path: Optional SLOW5/BLOW5 file path
#         max_reads: Optional maximum number of reads to process

#     Yields:
#         Tuple of (read_id, alignments) for each read

#     Example:
#         >>> for read_id, alignments in eventalign("reads.bam", "ref.fa", "signals.blow5"):
#         ...     print(f"{read_id}: {len(alignments)} events")
#     """
#     with EventAligner(bam_path, fasta_path, slow5_path) as aligner:
#         yield from aligner.align_reads(max_reads=max_reads)


# __all__ = [
#     "EventAligner",
#     "eventalign",
#     "_F5C_AVAILABLE",
# ]


# if __name__ == "__main__":
#     # Simple test
#     if _F5C_AVAILABLE:
#         print("✓ f5c module is available")
#         print("  Use EventAligner class to perform event-to-kmer alignments")
#     else:
#         print("✗ f5c module not available")
#         print("  Install with: pip install -e .")
