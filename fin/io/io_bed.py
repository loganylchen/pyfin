"""
BED file format parser
Supports BED3, BED6, BED12 formats
"""

from typing import List, Optional, Tuple, Generator
from pathlib import Path
import logging


logger = logging.getLogger(__name__)


class BEDFeature:
    """Represents a feature/interval in a BED file"""

    def __init__(self, chrom: str, start: int, end: int, name: Optional[str] = None,
                 score: Optional[float] = None, strand: Optional[str] = None,
                 thick_start: Optional[int] = None, thick_end: Optional[int] = None,
                 item_rgb: Optional[str] = None, block_count: Optional[int] = None,
                 block_sizes: Optional[List[int]] = None, block_starts: Optional[List[int]] = None,
                 extra_fields: Optional[List[str]] = None):
        """
        Initialize BED feature

        Args:
            chrom: Chromosome name
            start: Start position (0-based)
            end: End position (exclusive)
            name: Feature name
            score: Score (0-1000)
            strand: Strand (+, -, or .)
            thick_start: Thick start position
            thick_end: Thick end position
            item_rgb: RGB color values (r,g,b)
            block_count: Number of blocks/exons
            block_sizes: Block sizes
            block_starts: Block start positions
            extra_fields: Extra fields for custom BED formats
        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.thick_start = thick_start
        self.thick_end = thick_end
        self.item_rgb = item_rgb
        self.block_count = block_count
        self.block_sizes = block_sizes or []
        self.block_starts = block_starts or []
        self.extra_fields = extra_fields or []

    @property
    def length(self) -> int:
        """Get feature length"""
        return self.end - self.start

    @property
    def is_bed12(self) -> bool:
        """Check if this is a BED12 feature"""
        return self.block_count is not None and self.block_count > 0

    def __str__(self) -> str:
        """String representation (BED format)"""
        fields = [self.chrom, str(self.start), str(self.end)]

        if self.name is not None:
            fields.append(self.name)
        if self.score is not None:
            fields.append(str(self.score))
        if self.strand is not None:
            fields.append(self.strand)
        if self.thick_start is not None:
            fields.append(str(self.thick_start))
        if self.thick_end is not None:
            fields.append(str(self.thick_end))
        if self.item_rgb is not None:
            fields.append(self.item_rgb)
        if self.block_count is not None:
            fields.append(str(self.block_count))
        if self.block_sizes:
            fields.append(','.join(map(str, self.block_sizes)) + ',' if self.block_sizes else '0,')
        if self.block_starts:
            fields.append(','.join(map(str, self.block_starts)) + ',' if self.block_starts else '0,')
        if self.extra_fields:
            fields.extend(self.extra_fields)

        return '\t'.join(fields)


class BEDReader:
    """Reader for BED files"""

    def __init__(self, file_path: str):
        """
        Initialize BED reader

        Args:
            file_path: Path to the BED file
        """
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"BED file not found: {file_path}")

        self._file = None
        self._feature_count = None

    def __enter__(self):
        """Context manager entry"""
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.close()

    def open(self):
        """Open the BED file"""
        try:
            self._file = open(self.file_path, 'r')
            logger.info(f"Opened BED file: {self.file_path}")
        except Exception as e:
            logger.error(f"Failed to open BED file {self.file_path}: {e}")
            raise

    def close(self):
        """Close the BED file"""
        if self._file is not None:
            self._file.close()
            self._file = None
            logger.info(f"Closed BED file: {self.file_path}")

    def _parse_line(self, line: str) -> Optional[BEDFeature]:
        """
        Parse a single BED line

        Args:
            line: Line from BED file

        Returns:
            BEDFeature object or None if invalid
        """
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('browser ') or line.startswith('track '):
            return None

        fields = line.split('\t')
        if len(fields) < 3:
            logger.warning(f"Invalid BED line (minimum 3 fields required): {line}")
            return None

        try:
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            # Initialize optional fields
            name = fields[3] if len(fields) > 3 else None
            score = float(fields[4]) if len(fields) > 4 and fields[4] != '.' else None
            strand = fields[5] if len(fields) > 5 and fields[5] in ['+', '-', '.'] else None
            thick_start = int(fields[6]) if len(fields) > 6 else None
            thick_end = int(fields[7]) if len(fields) > 7 else None
            item_rgb = fields[8] if len(fields) > 8 else None
            block_count = int(fields[9]) if len(fields) > 9 else None

            block_sizes = None
            if len(fields) > 10 and fields[10] and fields[10] != '.':
                block_sizes = [int(x) for x in fields[10].rstrip(',').split(',') if x]

            block_starts = None
            if len(fields) > 11 and fields[11] and fields[11] != '.':
                block_starts = [int(x) for x in fields[11].rstrip(',').split(',') if x]

            extra_fields = fields[12:] if len(fields) > 12 else None

            return BEDFeature(
                chrom=chrom, start=start, end=end, name=name,
                score=score, strand=strand, thick_start=thick_start,
                thick_end=thick_end, item_rgb=item_rgb, block_count=block_count,
                block_sizes=block_sizes, block_starts=block_starts, extra_fields=extra_fields
            )

        except (ValueError, IndexError) as e:
            logger.warning(f"Failed to parse BED line: {line}\nError: {e}")
            return None

    def get_features(self, chrom: Optional[str] = None,
                     start: Optional[int] = None,
                     end: Optional[int] = None,
                     strand: Optional[str] = None) -> List[BEDFeature]:
        """
        Get all features, optionally filtered by region

        Args:
            chrom: Chromosome to filter by
            start: Start position to filter by
            end: End position to filter by
            strand: Strand to filter by (+, -, or .)

        Returns:
            List of BEDFeature objects
        """
        if self._file is None:
            raise RuntimeError("BED file not opened. Call open() first.")

        # Reset file pointer to beginning
        self._file.seek(0)

        features = []
        for line_num, line in enumerate(self._file, 1):
            try:
                feature = self._parse_line(line)
                if feature is None:
                    continue

                # Apply filters
                if chrom and feature.chrom != chrom:
                    continue
                if start is not None and feature.end <= start:
                    continue
                if end is not None and feature.start >= end:
                    continue
                if strand and feature.strand != strand:
                    continue

                features.append(feature)

            except Exception as e:
                logger.warning(f"Error parsing line {line_num}: {e}")

        return features

    def iterate_features(self) -> Generator[BEDFeature, None, None]:
        """
        Generator to iterate through all features

        Yields:
            BEDFeature objects
        """
        if self._file is None:
            raise RuntimeError("BED file not opened. Call open() first.")

        # Reset file pointer to beginning
        self._file.seek(0)

        for line_num, line in enumerate(self._file, 1):
            try:
                feature = self._parse_line(line)
                if feature is not None:
                    yield feature
            except Exception as e:
                logger.warning(f"Error parsing line {line_num}: {e}")

    def get_feature_count(self) -> int:
        """
        Get total number of features

        Returns:
            Number of features
        """
        if self._feature_count is None:
            if self._file is None:
                raise RuntimeError("BED file not opened. Call open() first.")

            self._file.seek(0)
            count = 0
            for line in self._file:
                if self._parse_line(line) is not None:
                    count += 1
            self._feature_count = count

        return self._feature_count

    def get_chromosomes(self) -> List[str]:
        """
        Get list of chromosomes in the file

        Returns:
            List of chromosome names
        """
        chroms = set()
        for feature in self.iterate_features():
            chroms.add(feature.chrom)
        return sorted(list(chroms))

    def get_strand_specific_features(self) -> Tuple[List[BEDFeature], List[BEDFeature], List[BEDFeature]]:
        """
        Get features separated by strand

        Returns:
            Tuple of (plus_strand_features, minus_strand_features, unstranded_features)
        """
        plus, minus, unstranded = [], [], []

        for feature in self.iterate_features():
            if feature.strand == '+':
                plus.append(feature)
            elif feature.strand == '-':
                minus.append(feature)
            else:
                unstranded.append(feature)

        return plus, minus, unstranded

    def overlaps(self, feature1: BEDFeature, feature2: BEDFeature) -> bool:
        """
        Check if two features overlap

        Args:
            feature1: First feature
            feature2: Second feature

        Returns:
            True if features overlap
        """
        return (feature1.chrom == feature2.chrom and
                feature1.start < feature2.end and
                feature2.start < feature1.end)

    def merge_features(self, features: List[BEDFeature]) -> List[BEDFeature]:
        """
        Merge overlapping features

        Args:
            features: List of features to merge

        Returns:
            List of merged features
        """
        if not features:
            return []

        # Sort by chromosome and start position
        features = sorted(features, key=lambda f: (f.chrom, f.start))

        merged = []
        current = features[0]

        for feature in features[1:]:
            if (feature.chrom == current.chrom and
                feature.start <= current.end):  # Overlapping
                # Merge: extend end if needed
                current = BEDFeature(
                    chrom=current.chrom,
                    start=current.start,
                    end=max(current.end, feature.end),
                    name=current.name,
                    strand=current.strand
                )
            else:
                merged.append(current)
                current = feature

        merged.append(current)
        return merged

    @staticmethod
    def write_bed(features: List[BEDFeature], output_path: str, header: Optional[str] = None):
        """
        Write features to a BED file

        Args:
            features: List of BEDFeature objects
            output_path: Output file path
            header: Optional header line
        """
        try:
            with open(output_path, 'w') as f:
                if header:
                    if not header.startswith('track ') and not header.startswith('browser '):
                        f.write(f"# {header}\n")
                    else:
                        f.write(f"{header}\n")

                for feature in features:
                    f.write(str(feature) + '\n')

            logger.info(f"Wrote {len(features)} features to {output_path}")

        except Exception as e:
            logger.error(f"Failed to write BED file {output_path}: {e}")
            raise
