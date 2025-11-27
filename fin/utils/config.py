"""
Configuration management for FIN isoform detection pipeline.

Provides YAML-based configuration system for all algorithm parameters,
including eventalign settings, DTW parameters, and quality thresholds.
"""

import yaml
from dataclasses import dataclass, field
from typing import Optional, Dict, Any, Tuple
from pathlib import Path


@dataclass
class EventalignConfig:
    """Configuration for eventalign parameters."""
    min_event_duration: int = 1
    adaptive_bandwidth: bool = True
    min_eventalign_score: float = 0.8
    max_gap: int = 2
    min_events_per_kmer: int = 1

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "min_event_duration": self.min_event_duration,
            "adaptive_bandwidth": self.adaptive_bandwidth,
            "min_eventalign_score": self.min_eventalign_score,
            "max_gap": self.max_gap,
            "min_events_per_kmer": self.min_events_per_kmer
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'EventalignConfig':
        """Create from dictionary."""
        return cls(
            min_event_duration=data.get("min_event_duration", 1),
            adaptive_bandwidth=data.get("adaptive_bandwidth", True),
            min_eventalign_score=data.get("min_eventalign_score", 0.8),
            max_gap=data.get("max_gap", 2),
            min_events_per_kmer=data.get("min_events_per_kmer", 1)
        )


@dataclass
class DTWConfig:
    """Configuration for DTW (Dynamic Time Warping) parameters."""
    use_cuda: bool = True
    open_start: bool = False
    open_end: bool = False
    normalize_distance: bool = True
    batch_size: int = 1000
    max_distance: float = 10.0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "use_cuda": self.use_cuda,
            "open_start": self.open_start,
            "open_end": self.open_end,
            "normalize_distance": self.normalize_distance,
            "batch_size": self.batch_size,
            "max_distance": self.max_distance
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'DTWConfig':
        """Create from dictionary."""
        return cls(
            use_cuda=data.get("use_cuda", True),
            open_start=data.get("open_start", False),
            open_end=data.get("open_end", False),
            normalize_distance=data.get("normalize_distance", True),
            batch_size=data.get("batch_size", 1000),
            max_distance=data.get("max_distance", 10.0)
        )


@dataclass
class AlgorithmConfig:
    """Configuration for main algorithm parameters."""
    min_completeness_threshold: float = 0.80
    dtw_similarity_threshold: float = 0.15
    min_read_support: int = 5
    min_isoform_length: int = 200
    max_intron_size: int = 100000
    cluster_consistency_threshold: float = 0.6

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "min_completeness_threshold": self.min_completeness_threshold,
            "dtw_similarity_threshold": self.dtw_similarity_threshold,
            "min_read_support": self.min_read_support,
            "min_isoform_length": self.min_isoform_length,
            "max_intron_size": self.max_intron_size,
            "cluster_consistency_threshold": self.cluster_consistency_threshold
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'AlgorithmConfig':
        """Create from dictionary."""
        return cls(
            min_completeness_threshold=data.get("min_completeness_threshold", 0.80),
            dtw_similarity_threshold=data.get("dtw_similarity_threshold", 0.15),
            min_read_support=data.get("min_read_support", 5),
            min_isoform_length=data.get("min_isoform_length", 200),
            max_intron_size=data.get("max_intron_size", 100000),
            cluster_consistency_threshold=data.get("cluster_consistency_threshold", 0.6)
        )


@dataclass
class OutputConfig:
    """Configuration for output settings."""
    output_dir: str = "outputs/"
    write_intermediate: bool = True
    compress_output: bool = False
    bed12_format: bool = True
    add_metadata: bool = True
    precision: int = 3

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "output_dir": self.output_dir,
            "write_intermediate": self.write_intermediate,
            "compress_output": self.compress_output,
            "bed12_format": self.bed12_format,
            "add_metadata": self.add_metadata,
            "precision": self.precision
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'OutputConfig':
        """Create from dictionary."""
        return cls(
            output_dir=data.get("output_dir", "outputs/"),
            write_intermediate=data.get("write_intermediate", True),
            compress_output=data.get("compress_output", False),
            bed12_format=data.get("bed12_format", True),
            add_metadata=data.get("add_metadata", True),
            precision=data.get("precision", 3)
        )


@dataclass
class PipelineConfig:
    """Main configuration container for the complete pipeline."""
    algorithm: AlgorithmConfig = field(default_factory=AlgorithmConfig)
    eventalign: EventalignConfig = field(default_factory=EventalignConfig)
    dtw: DTWConfig = field(default_factory=DTWConfig)
    output: OutputConfig = field(default_factory=OutputConfig)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to nested dictionary."""
        return {
            "algorithm": self.algorithm.to_dict(),
            "eventalign": self.eventalign.to_dict(),
            "dtw": self.dtw.to_dict(),
            "output": self.output.to_dict()
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'PipelineConfig':
        """Create from nested dictionary."""
        return cls(
            algorithm=AlgorithmConfig.from_dict(data.get("algorithm", {})),
            eventalign=EventalignConfig.from_dict(data.get("eventalign", {})),
            dtw=DTWConfig.from_dict(data.get("dtw", {})),
            output=OutputConfig.from_dict(data.get("output", {}))
        )

    def to_yaml(self, output_path: str):
        """Write configuration to YAML file."""
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False, sort_keys=False)
        return output_path

    @classmethod
    def from_yaml(cls, yaml_path: str) -> 'PipelineConfig':
        """Load configuration from YAML file."""
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)
        return cls.from_dict(data or {})


class ConfigManager:
    """Manager class for configuration operations."""

    @staticmethod
    def create_default_config(output_path: Optional[str] = None) -> PipelineConfig:
        """
        Create configuration with default values.

        Args:
            output_path: If provided, write to this path

        Returns:
            PipelineConfig object
        """
        config = PipelineConfig()

        if output_path:
            config.to_yaml(output_path)

        return config

    @staticmethod
    def load_config(config_path: str) -> PipelineConfig:
        """
        Load configuration from file.

        Args:
            config_path: Path to YAML configuration file

        Returns:
            PipelineConfig object

        Raises:
            FileNotFoundError: If config file doesn't exist
        """
        return PipelineConfig.from_yaml(config_path)

    @staticmethod
    def validate_config(config: PipelineConfig) -> Tuple[bool, str]:
        """
        Validate configuration parameters.

        Args:
            config: PipelineConfig object to validate

        Returns:
            Tuple of (is_valid, error_message)
        """
        # Check completeness threshold
        if not (0.0 <= config.algorithm.min_completeness_threshold <= 1.0):
            return False, "min_completeness_threshold must be between 0 and 1"

        # Check DTW threshold
        if not (0.0 <= config.algorithm.dtw_similarity_threshold <= 1.0):
            return False, "dtw_similarity_threshold must be between 0 and 1"

        # Check read support
        if config.algorithm.min_read_support < 1:
            return False, "min_read_support must be at least 1"

        # Check isoform length
        if config.algorithm.min_isoform_length < 1:
            return False, "min_isoform_length must be positive"

        # Check eventalign parameters
        if config.eventalign.min_event_duration < 1:
            return False, "min_event_duration must be at least 1"

        return True, "Configuration is valid"

    @staticmethod
    def merge_configs(base_config: PipelineConfig, override_dict: Dict[str, Any]) -> PipelineConfig:
        """
        Merge base configuration with override dictionary.

        Args:
            base_config: Base configuration
            override_dict: Nested dictionary with override values

        Returns:
            New PipelineConfig with merged values
        """
        base_dict = base_config.to_dict()

        # Merge nested dictionaries
        def merge_dicts(d1, d2):
            for key, value in d2.items():
                if key in d1 and isinstance(d1[key], dict) and isinstance(value, dict):
                    merge_dicts(d1[key], value)
                else:
                    d1[key] = value

        merge_dicts(base_dict, override_dict)
        return PipelineConfig.from_dict(base_dict)


# Example configuration file template
DEFAULT_CONFIG_TEMPLATE = """
# Algorithm parameters
algorithm:
  min_completeness_threshold: 0.80   # Minimum eventalign completeness score (0-1)
  dtw_similarity_threshold: 0.15     # Maximum normalized DTW distance (0-1)
  min_read_support: 5                  # Minimum reads to validate isoform
  min_isoform_length: 200              # Minimum transcript length (bp)
  max_intron_size: 100000              # Maximum allowed intron size (bp)
  cluster_consistency_threshold: 0.6   # Minimum cluster consistency (0-1)

# Eventalign parameters
eventalign:
  min_event_duration: 1                # Minimum event duration (samples)
  adaptive_bandwidth: true             # Use adaptive banding
  min_eventalign_score: 0.8            # Minimum alignment score (0-1)
  max_gap: 2                          # Maximum gap between events
  min_events_per_kmer: 1              # Minimum events per k-mer

# DTW parameters
dtw:
  use_cuda: true                       # Enable GPU acceleration
  open_start: false                    # Allow open start alignment
  open_end: false                      # Allow open end alignment
  normalize_distance: true            # Normalize by sequence length
  batch_size: 1000                   # Batch size for DTW computation
  max_distance: 10.0                  # Maximum allowed DTW distance

# Output settings
output:
  output_dir: "./outputs"              # Output directory
  write_intermediate: true            # Write intermediate results
  compress_output: false              # Compress output files
  bed12_format: true                  # Output read assignments in BED12
  add_metadata: true                  # Add metadata to output files
  precision: 3                        # Decimal precision for scores
"""


if __name__ == "__main__":
    # Create example config
    import tempfile

    config = ConfigManager.create_default_config()

    # Print as YAML
    print("# Example configuration:")
    print(DEFAULT_CONFIG_TEMPLATE)

    # Validate
    is_valid, message = ConfigManager.validate_config(config)
    print(f"\nConfiguration validation: {message}")
