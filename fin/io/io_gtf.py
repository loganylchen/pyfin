"""
GTF/GFF file format parser and writer
"""

from __future__ import annotations

from typing import List, Dict, Optional, Any, Generator, TYPE_CHECKING
from pathlib import Path
import logging

if TYPE_CHECKING:
    from fin.analysis.quantification import QuantResult


logger = logging.getLogger(__name__)


class GTFTranscript:
    """Represents a transcript with exons and optional CDS features"""

    def __init__(self, transcript_id: str, gene_id: str, chrom: str,
                 strand: str, source: str = 'unknown', score: Optional[float] = None,
                 frame: Optional[int] = None, attributes: Optional[Dict[str, str]] = None):
        """
        Initialize transcript

        Args:
            transcript_id: Transcript ID
            gene_id: Gene ID
            chrom: Chromosome name
            strand: Strand (+ or -)
            source: Source/annotation origin
            score: Score value
            frame: Frame information (0, 1, 2, or .)
            attributes: Additional attributes
        """
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.chrom = chrom
        self.strand = strand
        self.source = source
        self.score = score
        self.frame = frame
        self.attributes = attributes or {}

        self.exons = []  # List of (start, end) tuples
        self.cds = []    # List of (start, end) tuples
        self.start_codon = None
        self.stop_codon = None
        self.utr5 = []   # 5' UTR regions
        self.utr3 = []   # 3' UTR regions

    @property
    def start(self) -> int:
        """Get transcript start (minimum position of all features)"""
        all_starts = []
        if self.exons:
            all_starts.extend([exon[0] for exon in self.exons])
        if self.cds:
            all_starts.extend([cds[0] for cds in self.cds])
        return min(all_starts) if all_starts else 0

    @property
    def end(self) -> int:
        """Get transcript end (maximum position of all features)"""
        all_ends = []
        if self.exons:
            all_ends.extend([exon[1] for exon in self.exons])
        if self.cds:
            all_ends.extend([cds[1] for cds in self.cds])
        return max(all_ends) if all_ends else 0

    @property
    def length(self) -> int:
        """Get total transcript length including introns"""
        return self.end - self.start

    @property
    def cds_length(self) -> int:
        """Get total CDS length"""
        return sum(end - start for start, end in self.cds)

    @property
    def exon_length(self) -> int:
        """Get total exon length"""
        return sum(end - start for start, end in self.exons)

    @property
    def num_exons(self) -> int:
        """Get number of exons"""
        return len(self.exons)

    def sort_features(self):
        """Sort exons and CDS by position"""
        self.exons.sort()
        self.cds.sort()

    def add_exon(self, start: int, end: int):
        """Add an exon"""
        self.exons.append((start, end))

    def add_cds(self, start: int, end: int):
        """Add a CDS region"""
        self.cds.append((start, end))

    def get_spliced_sequence(self, fasta_seq: str) -> str:
        """
        Get spliced sequence for this transcript

        Args:
            fasta_seq: Full chromosome sequence as string

        Returns:
            Spliced transcript sequence
        """
        if not self.exons:
            return ""

        sequence = ""
        for start, end in self.exons:
            sequence += fasta_seq[start:end]

        # Reverse complement for negative strand
        if self.strand == '-':
            sequence = self._reverse_complement(sequence)

        return sequence

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        """Get reverse complement of a sequence"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                     'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
        return "".join(complement.get(base, 'N') for base in reversed(seq))

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'transcript_id': self.transcript_id,
            'gene_id': self.gene_id,
            'chrom': self.chrom,
            'strand': self.strand,
            'start': self.start,
            'end': self.end,
            'length': self.length,
            'source': self.source,
            'score': self.score,
            'frame': self.frame,
            'attributes': self.attributes,
            'exons': self.exons,
            'num_exons': self.num_exons,
            'exon_length': self.exon_length,
            'cds': self.cds,
            'cds_length': self.cds_length,
        }


class GTFGene:
    """Represents a gene with multiple transcripts"""

    def __init__(self, gene_id: str, chrom: str, strand: str, source: str = 'unknown',
                 score: Optional[float] = None, attributes: Optional[Dict[str, str]] = None):
        """
        Initialize gene

        Args:
            gene_id: Gene ID
            chrom: Chromosome name
            strand: Strand (+ or -)
            source: Source/annotation origin
            score: Score value
            attributes: Additional attributes
        """
        self.gene_id = gene_id
        self.chrom = chrom
        self.strand = strand
        self.source = source
        self.score = score
        self.attributes = attributes or {}
        self.transcripts = {}  # transcript_id -> GTFTranscript

    @property
    def start(self) -> int:
        """Get gene start (minimum position of all transcripts)"""
        if not self.transcripts:
            return 0
        return min(tx.start for tx in self.transcripts.values())

    @property
    def end(self) -> int:
        """Get gene end (maximum position of all transcripts)"""
        if not self.transcripts:
            return 0
        return max(tx.end for tx in self.transcripts.values())

    @property
    def length(self) -> int:
        """Get gene length"""
        return self.end - self.start

    @property
    def num_transcripts(self) -> int:
        """Get number of transcripts"""
        return len(self.transcripts)

    def add_transcript(self, transcript: GTFTranscript):
        """Add a transcript to the gene"""
        self.transcripts[transcript.transcript_id] = transcript

    def get_cds_length(self) -> int:
        """Get max CDS length among all transcripts"""
        if not self.transcripts:
            return 0
        return max(tx.cds_length for tx in self.transcripts.values())

    def get_exon_numbers(self) -> List[int]:
        """Get list of exon numbers for all transcripts"""
        return [tx.num_exons for tx in self.transcripts.values()]

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'gene_id': self.gene_id,
            'chrom': self.chrom,
            'strand': self.strand,
            'start': self.start,
            'end': self.end,
            'length': self.length,
            'source': self.source,
            'score': self.score,
            'attributes': self.attributes,
            'num_transcripts': self.num_transcripts,
            'transcript_ids': list(self.transcripts.keys()),
        }


class GTFReader:
    """Reader for GTF/GFF files"""

    def __init__(self, file_path: str):
        """
        Initialize GTF/GFF reader

        Args:
            file_path: Path to GTF/GFF file
        """
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"GTF/GFF file not found: {file_path}")

        self._file = None
        self.genes = {}  # gene_id -> GTFGene
        self.transcripts = {}  # transcript_id -> GTFTranscript
        self._parsed = False

    def __enter__(self):
        """Context manager entry"""
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.close()

    def open(self):
        """Open GTF/GFF file"""
        try:
            self._file = open(self.file_path, 'r')
            logger.info(f"Opened GTF/GFF file: {self.file_path}")
        except Exception as e:
            logger.error(f"Failed to open GTF/GFF file {self.file_path}: {e}")
            raise

    def close(self):
        """Close file"""
        if self._file is not None:
            self._file.close()
            self._file = None
            logger.info(f"Closed GTF/GFF file: {self.file_path}")

    def _parse_attributes(self, attr_str: str) -> Dict[str, str]:
        """
        Parse GTF/GFF attribute field

        Args:
            attr_str: Attribute string (key "value"; key2 "value2"; ...)

        Returns:
            Dictionary of attributes
        """
        if not attr_str or attr_str == '.':
            return {}

        attrs = {}

        # Handle different formats (GTF with space separation, GFF with =)
        if '; ' in attr_str or (attr_str.count(';') > 1 and '"' in attr_str):
            # GTF format: key "value"; key2 "value2";
            for attr in attr_str.rstrip(';').split('; '):
                if not attr.strip():
                    continue
                parts = attr.strip().split(' ', 1)
                if len(parts) == 2:
                    key = parts[0].strip()
                    value = parts[1].strip().strip('"')
                    attrs[key] = value
        else:
            # GFF3 format: key=value; key2=value2
            for attr in attr_str.rstrip(';').split(';'):
                if not attr.strip():
                    continue
                if '=' in attr:
                    key, value = attr.strip().split('=', 1)
                    attrs[key.strip()] = value.strip().strip('"')
                elif ' ' in attr:
                    # Fallback to space-separated
                    parts = attr.strip().split(' ', 1)
                    if len(parts) == 2:
                        key = parts[0].strip()
                        value = parts[1].strip().strip('"')
                        attrs[key] = value

        return attrs

    def _extract_id(self, attributes: Dict[str, str], id_keys: List[str]) -> Optional[str]:
        """
        Extract ID from attributes using multiple possible keys

        Args:
            attributes: Attribute dictionary
            id_keys: List of possible ID keys to try

        Returns:
            ID string or None if not found
        """
        for key in id_keys:
            if key in attributes:
                return attributes[key]
        return None

    def parse(self):
        """
        Parse entire GTF/GFF file and build gene/transcript structure
        """
        if self._file is None:
            raise RuntimeError("GTF/GFF file not opened. Call open() first.")

        self._file.seek(0)

        # Temporary storage for transcript components
        transcript_parts = {}  # transcript_id -> {'exons': [], 'cds': [], ...}

        # Common attribute keys to look for
        gene_id_keys = ['gene_id', 'geneid', 'geneID', 'ID']
        transcript_id_keys = ['transcript_id', 'transcriptid', 'transcriptID', 'Parent']

        # Keep track of current gene/transcript
        current_gene_id = None
        current_transcript_id = None

        for line_num, line in enumerate(self._file, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            try:
                if line.startswith('#'):
                    logger.warning('Comment line, skipped!')
                # Parse GTF/GFF line
                fields = line.split('\t')
                if len(fields) < 9:
                    logger.warning(f"Invalid GTF/GFF line at {line_num}: not enough fields")
                    continue

                chrom, source, feature_type, start_str, end_str, score_str, strand, frame_str, attr_str = fields[:9]

                # Convert positions (GTF is 1-based, convert to 0-based)
                start = int(start_str) - 1
                end = int(end_str)
                score = float(score_str) if score_str != '.' else None
                frame = int(frame_str) if frame_str.isdigit() else None

                # Parse attributes
                attributes = self._parse_attributes(attr_str)

                # Extract IDs
                gene_id = self._extract_id(attributes, gene_id_keys)
                transcript_id = self._extract_id(attributes, transcript_id_keys)

                if not gene_id and not transcript_id:
                    logger.warning(f"No gene or transcript ID found at line {line_num}")
                    continue

                # Use gene_id as transcript_id if only gene_id exists (simpler annotation)
                if gene_id and not transcript_id:
                    transcript_id = gene_id
                elif transcript_id and not gene_id:
                    # Try to find gene_id in attributes or use transcript_id as gene_id
                    for key in ['gene_id', 'geneid', 'geneID', 'gene', 'ID']:
                        gene_id = attributes.get(key, transcript_id)
                        if gene_id:
                            break

                # Create/update gene
                if gene_id not in self.genes:
                    self.genes[gene_id] = GTFGene(
                        gene_id=gene_id, chrom=chrom, strand=strand,
                        source=source, score=score, attributes=attributes
                    )

                # Create/update transcript
                if transcript_id not in self.transcripts:
                    self.transcripts[transcript_id] = GTFTranscript(
                        transcript_id=transcript_id, gene_id=gene_id, chrom=chrom,
                        strand=strand, source=source, score=score, frame=frame,
                        attributes=attributes
                    )
                    self.genes[gene_id].add_transcript(self.transcripts[transcript_id])

                transcript_obj = self.transcripts[transcript_id]

                # Add features based on feature_type
                feature_type_lower = feature_type.lower()

                if 'exon' in feature_type_lower:
                    transcript_obj.add_exon(start, end)
                elif 'cds' in feature_type_lower:
                    transcript_obj.add_cds(start, end)
                elif 'utr' in feature_type_lower:
                    if '5' in feature_type_lower or 'five' in feature_type_lower:
                        transcript_obj.utr5.append((start, end))
                    elif '3' in feature_type_lower or 'three' in feature_type_lower:
                        transcript_obj.utr3.append((start, end))
                elif 'start_codon' in feature_type_lower:
                    transcript_obj.start_codon = (start, end)
                elif 'stop_codon' in feature_type_lower:
                    transcript_obj.stop_codon = (start, end)

            except Exception as e:
                logger.warning(f"Error parsing line {line_num}: {e}: {fields[:9]}")
                continue

        # Sort features for all transcripts
        for transcript in self.transcripts.values():
            transcript.sort_features()

        self._parsed = True
        logger.info(f"Parsed {len(self.genes)} genes and {len(self.transcripts)} transcripts")

    def get_gene(self, gene_id: str) -> Optional[GTFGene]:
        """Get gene by ID"""
        return self.genes.get(gene_id)

    def get_transcript(self, transcript_id: str) -> Optional[GTFTranscript]:
        """Get transcript by ID"""
        return self.transcripts.get(transcript_id)

    def get_genes_on_chromosome(self, chrom: str) -> List[GTFGene]:
        """
        Get all genes on a specific chromosome

        Args:
            chrom: Chromosome name

        Returns:
            List of GTFGene objects
        """
        return [gene for gene in self.genes.values() if gene.chrom == chrom]

    def get_transcripts_on_chromosome(self, chrom: str) -> List[GTFTranscript]:
        """
        Get all transcripts on a specific chromosome

        Args:
            chrom: Chromosome name

        Returns:
            List of GTFTranscript objects
        """
        return [tx for tx in self.transcripts.values() if tx.chrom == chrom]

    def get_genes_in_region(self, chrom: str, start: int, end: int) -> List[GTFGene]:
        """
        Get genes overlapping a region

        Args:
            chrom: Chromosome name
            start: Start position (0-based)
            end: End position

        Returns:
            List of GTFGene objects
        """
        genes_in_region = []
        for gene in self.genes.values():
            if gene.chrom == chrom and gene.start < end and gene.end > start:
                genes_in_region.append(gene)
        return genes_in_region

    def get_transcripts_in_region(self, chrom: str, start: int, end: int) -> List[GTFTranscript]:
        """
        Get transcripts overlapping a region

        Args:
            chrom: Chromosome name
            start: Start position (0-based)
            end: End position

        Returns:
            List of GTFTranscript objects
        """
        txs_in_region = []
        for tx in self.transcripts.values():
            if tx.chrom == chrom and tx.start < end and tx.end > start:
                txs_in_region.append(tx)
        return txs_in_region

    def iterate_genes(self) -> Generator[GTFGene, None, None]:
        """Iterate over all genes"""
        for gene in self.genes.values():
            yield gene

    def iterate_transcripts(self) -> Generator[GTFTranscript, None, None]:
        """Iterate over all transcripts"""
        for transcript in self.transcripts.values():
            yield transcript

    def get_gene_stats(self) -> Dict[str, Any]:
        """
        Get statistics about genes in the file

        Returns:
            Dictionary with statistics
        """
        if not self.genes:
            return {}

        gene_lengths = [gene.length for gene in self.genes.values()]
        transcript_counts = [gene.num_transcripts for gene in self.genes.values()]

        chromosomes = set(gene.chrom for gene in self.genes.values())
        strands = set(gene.strand for gene in self.genes.values())

        return {
            'num_genes': len(self.genes),
            'num_transcripts': len(self.transcripts),
            'num_chromosomes': len(chromosomes),
            'strands': list(strands),
            'avg_gene_length': sum(gene_lengths) / len(gene_lengths),
            'median_gene_length': sorted(gene_lengths)[len(gene_lengths) // 2],
            'avg_transcripts_per_gene': sum(transcript_counts) / len(transcript_counts),
            'genes_per_chromosome': {chrom: sum(1 for g in self.genes.values() if g.chrom == chrom)
                                    for chrom in chromosomes},
        }


def write_gtf(
    results: Dict[str, "QuantResult"],
    output_path: str,
    source: str = "pyfin",
) -> None:
    """Write quantification results as a GTF file with gene/transcript/exon hierarchy.

    Args:
        results: Dict mapping candidate_id -> QuantResult (with genomic fields populated).
        output_path: Path to write the GTF file.
        source: Value for the GTF source column.
    """

    # Group transcripts by gene_id
    genes: Dict[str, List[QuantResult]] = {}
    for qr in results.values():
        gid = qr.gene_id or qr.candidate_id
        genes.setdefault(gid, []).append(qr)

    # Collect all records, sort by (chrom, start)
    sorted_genes = sorted(
        genes.items(),
        key=lambda item: (item[1][0].chrom, min(qr.start for qr in item[1])),
    )

    lines: List[str] = []
    for gene_id, transcripts in sorted_genes:
        transcripts.sort(key=lambda qr: qr.start)
        chrom = transcripts[0].chrom
        strand = transcripts[0].strand
        gene_start = min(qr.start for qr in transcripts)
        gene_end = max(qr.end for qr in transcripts)

        # Gene line (0-based internal -> 1-based GTF)
        lines.append(_gtf_line(
            chrom, source, "gene", gene_start + 1, gene_end, ".", strand, ".",
            f'gene_id "{gene_id}";',
        ))

        for qr in transcripts:
            tid = qr.candidate_id

            # Fusion chrom splitting: seqname uses chromA; add fusion_partner_chrom attr
            qr_chrom = qr.chrom
            fusion_attr = ""
            if getattr(qr, "source", None) == "fusion" and "::" in qr_chrom:
                chrom_a, chrom_b = qr_chrom.split("::", 1)
                qr_chrom = chrom_a
                fusion_attr = f'fusion_partner_chrom "{chrom_b}"; '

            coherence = getattr(qr, "coherence_score", 0.0)
            discrimination = getattr(qr, "discrimination_score", 0.0)
            combined = getattr(qr, "combined_score", 0.0)

            attrs = (
                f'gene_id "{gene_id}"; '
                f'transcript_id "{tid}"; '
                f'abundance "{qr.abundance:.4f}"; '
                f'confidence "{qr.confidence:.4f}"; '
                f'num_reads "{qr.num_assigned_reads}"; '
                f'transcript_source "{qr.source}"; '
                f'{fusion_attr}'
                f'coherence_score "{coherence:.4f}"; '
                f'discrimination_score "{discrimination:.4f}"; '
                f'combined_score "{combined:.4f}";'
            )
            lines.append(_gtf_line(
                qr_chrom, source, "transcript", qr.start + 1, qr.end,
                ".", strand, ".", attrs,
            ))

            for exon_num, (es, ee) in enumerate(qr.exons, 1):
                exon_attrs = (
                    f'gene_id "{gene_id}"; '
                    f'transcript_id "{tid}"; '
                    f'exon_number "{exon_num}";'
                )
                lines.append(_gtf_line(
                    qr_chrom, source, "exon", es + 1, ee,
                    ".", strand, ".", exon_attrs,
                ))

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        for line in lines:
            fh.write(line)
            fh.write("\n")

    logger.info("Wrote %d GTF lines to %s", len(lines), output_path)


def _gtf_line(
    seqname: str, source: str, feature: str,
    start: int, end: int, score: str, strand: str, frame: str,
    attributes: str,
) -> str:
    """Format a single GTF line."""
    return f"{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attributes}"
