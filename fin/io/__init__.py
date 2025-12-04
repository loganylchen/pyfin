"""
I/O module for handling various data formats in nanopore Direct RNA-seq analysis
"""

# Import core readers that should always be available
try:
    from .io_fasta import FASTAReader, FASTARecord
    _has_fasta = True
except ImportError:
    _has_fasta = False

# Import interval manager functions
from .interval_manager import (
    GenomicInterval,
    generate_isolated_intervals,
    extract_reads_for_interval,
    extract_annotation_for_interval,
    is_fusion_read
)

# Import optional readers with lazy loading
class _LazyImporter:
    """Lazy loader for optional dependencies"""

    def __init__(self, module_name, class_name, import_error_msg):
        self.module_name = module_name
        self.class_name = class_name
        self.import_error_msg = import_error_msg
        self._module = None

    def __call__(self, *args, **kwargs):
        if self._module is None:
            try:
                module = __import__(f"fin.io.{self.module_name}", fromlist=[self.class_name])
                self._module = getattr(module, self.class_name)
            except ImportError as e:
                raise ImportError(self.import_error_msg) from e
        return self._module(*args, **kwargs)

    def __getattr__(self, name):
        if self._module is None:
            try:
                module = __import__(f"fin.io.{self.module_name}", fromlist=[self.class_name])
                self._module = getattr(module, self.class_name)
            except ImportError as e:
                raise ImportError(self.import_error_msg) from e
        return getattr(self._module, name)

# Create lazy loaders
Fast5Reader = _LazyImporter("io_fast5", "Fast5Reader",
    "ont-fast5-api is required for Fast5 support. Install it with: pip install ont-fast5-api")

Slow5Reader = _LazyImporter("io_slow5", "Slow5Reader",
    "pyslow5 is required for Slow5 support. Install it with: pip install pyslow5")

Pod5Reader = _LazyImporter("io_pod5", "Pod5Reader",
    "pod5 is required for Pod5 support. Install it with: pip install pod5")

BamReader = _LazyImporter("io_bam", "BamReader",
    "pysam is required for BAM support. Install it with: pip install pysam")

GTFReader = _LazyImporter("io_gtf", "GTFReader",
    "No additional dependencies required for GTFReader")

BEDReader = _LazyImporter("io_bed", "BEDReader",
    "No additional dependencies required for BEDReader")

ReadSubsetManager = _LazyImporter("io_read_manager", "ReadSubsetManager",
    "pysam is required for ReadSubsetManager. Install it with: pip install pysam")

ReadSubset = _LazyImporter("io_read_manager", "ReadSubset",
    "No additional dependencies required for ReadSubset")

ReadDataBundle = _LazyImporter("io_read_manager", "ReadDataBundle",
    "No additional dependencies required for ReadDataBundle")

if _has_fasta:
    _create_subset_manager = _LazyImporter("io_read_manager", "create_subset_manager",
        "pysam is required for create_subset_manager. Install it with: pip install pysam")

def create_subset_manager(*args, **kwargs):
    """Create a subset manager (with proper import error handling)"""
    if not _has_fasta:
        raise ImportError("Cannot import create_subset_manager: FASTA support libraries not available")
    return _create_subset_manager(*args, **kwargs)

# Public API
__all__ = [
    "Fast5Reader",
    "Slow5Reader",
    "Pod5Reader",
    "BamReader",
    "GTFReader",
    "BEDReader",
    "GenomicInterval",
    "generate_isolated_intervals",
    "extract_reads_for_interval",
    "extract_annotation_for_interval",
    "is_fusion_read",
]

# Add classes that imported successfully
if _has_fasta:
    __all__.extend([
        "FASTAReader",
        "FASTARecord",
        "ReadSubsetManager",
        "ReadSubset",
        "ReadDataBundle",
        "create_subset_manager",
    ])
