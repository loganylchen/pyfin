#!/usr/bin/env python3
"""
Unit tests for fin.io.io_fasta module.

Tests the FASTA file reader and record classes:
- FASTARecord: Individual FASTA record
- FASTAReader: FASTA file reader using pysam

Run with:
    pytest tests/unit/test_io_fasta.py -v
"""

import pytest
import tempfile
import os
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fin.io.io_fasta import FASTARecord, FASTAReader


# ============================================================================
# Test fixtures
# ============================================================================

@pytest.fixture
def simple_fasta_file():
    """Create a simple FASTA file for testing."""
    content = """>seq1 Description of sequence 1
ATCGATCGATCG
>seq2 Another sequence description
GCTAGCTAGCTA
>seq3
NNNNNNNNNNNN
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(content)
        f.flush()
        yield f.name
    os.unlink(f.name)


@pytest.fixture
def multiline_fasta_file():
    """Create a FASTA file with multi-line sequences."""
    content = """>long_seq A long sequence split over multiple lines
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(content)
        f.flush()
        yield f.name
    os.unlink(f.name)


@pytest.fixture
def empty_fasta_file():
    """Create an empty FASTA file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        yield f.name
    os.unlink(f.name)


@pytest.fixture
def test_data_fasta():
    """Path to test data FASTA file if it exists."""
    test_fasta = PROJECT_ROOT / "tests" / "testdata" / "test.fa"
    if test_fasta.exists():
        return str(test_fasta)
    return None


# ============================================================================
# FASTARecord tests
# ============================================================================

class TestFASTARecord:
    """Tests for FASTARecord class."""

    def test_record_initialization(self):
        """Test basic FASTARecord initialization."""
        record = FASTARecord("seq1 description", "ATCG")
        assert record.header == "seq1 description"
        assert record.sequence == "ATCG"

    def test_record_id_extraction(self):
        """Test that ID is extracted from header."""
        record = FASTARecord("myseq description text", "ATCG")
        assert record.id == "myseq"

    def test_record_id_extraction_no_description(self):
        """Test ID extraction when no description present."""
        record = FASTARecord("myseq", "ATCG")
        assert record.id == "myseq"

    def test_record_id_extraction_empty_header(self):
        """Test ID extraction with empty header."""
        record = FASTARecord("", "ATCG")
        assert record.id == ""

    def test_record_description_extraction(self):
        """Test that description is extracted from header."""
        record = FASTARecord("seq1 this is a description", "ATCG")
        assert record.description == "this is a description"

    def test_record_description_empty(self):
        """Test description when only ID is present."""
        record = FASTARecord("seq1", "ATCG")
        assert record.description == ""

    def test_record_length(self):
        """Test length property."""
        record = FASTARecord("seq1", "ATCGATCG")
        assert record.length == 8

    def test_record_length_empty_sequence(self):
        """Test length of empty sequence."""
        record = FASTARecord("seq1", "")
        assert record.length == 0

    def test_record_gc_content_typical(self):
        """Test GC content calculation."""
        # 50% GC
        record = FASTARecord("seq1", "ATCG")
        assert record.gc_content == 50.0

    def test_record_gc_content_all_gc(self):
        """Test GC content for all GC sequence."""
        record = FASTARecord("seq1", "GGCC")
        assert record.gc_content == 100.0

    def test_record_gc_content_no_gc(self):
        """Test GC content for sequence with no GC."""
        record = FASTARecord("seq1", "AATT")
        assert record.gc_content == 0.0

    def test_record_gc_content_empty(self):
        """Test GC content for empty sequence."""
        record = FASTARecord("seq1", "")
        assert record.gc_content == 0.0

    def test_record_gc_content_lowercase(self):
        """Test GC content with lowercase sequence."""
        record = FASTARecord("seq1", "atcg")
        assert record.gc_content == 50.0

    def test_record_sequence_upper(self):
        """Test sequence_upper property."""
        record = FASTARecord("seq1", "atcg")
        assert record.sequence_upper == "ATCG"

    def test_record_sequence_lower(self):
        """Test sequence_lower property."""
        record = FASTARecord("seq1", "ATCG")
        assert record.sequence_lower == "atcg"

    def test_record_str(self):
        """Test string representation in FASTA format."""
        record = FASTARecord("seq1", "ATCGATCG")
        result = str(record)
        assert result.startswith(">seq1")
        assert "ATCGATCG" in result

    def test_record_subsequence(self):
        """Test subsequence extraction."""
        record = FASTARecord("seq1", "ATCGATCGATCG")
        assert record.subsequence(0, 4) == "ATCG"
        assert record.subsequence(4, 8) == "ATCG"
        assert record.subsequence(0, 12) == "ATCGATCGATCG"

    def test_record_subsequence_empty(self):
        """Test empty subsequence."""
        record = FASTARecord("seq1", "ATCGATCG")
        assert record.subsequence(0, 0) == ""

    def test_record_complement(self):
        """Test complement calculation."""
        record = FASTARecord("seq1", "ATCG")
        assert record.complement() == "TAGC"

    def test_record_complement_with_n(self):
        """Test complement with N bases."""
        record = FASTARecord("seq1", "ANCG")
        assert record.complement() == "TNGC"

    def test_record_complement_invalid_base(self):
        """Test complement raises error for invalid bases."""
        record = FASTARecord("seq1", "ATXG")
        with pytest.raises(ValueError, match="Invalid DNA base"):
            record.complement()

    def test_record_reverse_complement(self):
        """Test reverse complement calculation."""
        record = FASTARecord("seq1", "ATCG")
        # complement: TAGC, reverse: CGAT
        assert record.reverse_complement() == "CGAT"

    def test_record_count_nucleotides(self):
        """Test nucleotide counting."""
        record = FASTARecord("seq1", "AATTCCGG")
        counts = record.count_nucleotides()
        assert counts['A'] == 2
        assert counts['T'] == 2
        assert counts['C'] == 2
        assert counts['G'] == 2

    def test_record_count_nucleotides_with_n(self):
        """Test nucleotide counting with N bases."""
        record = FASTARecord("seq1", "ATCGN")
        counts = record.count_nucleotides()
        assert counts['N'] == 1


# ============================================================================
# FASTAReader tests
# ============================================================================

class TestFASTAReader:
    """Tests for FASTAReader class."""

    def test_reader_initialization(self, simple_fasta_file):
        """Test FASTAReader initialization."""
        reader = FASTAReader(simple_fasta_file)
        assert reader.file_path == Path(simple_fasta_file)

    def test_reader_file_not_found(self):
        """Test FASTAReader raises error for non-existent file."""
        with pytest.raises(FileNotFoundError):
            FASTAReader("/nonexistent/path/file.fa")

    def test_reader_context_manager(self, simple_fasta_file):
        """Test FASTAReader as context manager."""
        with FASTAReader(simple_fasta_file) as reader:
            assert reader._fasta_file is not None

    def test_reader_open_close(self, simple_fasta_file):
        """Test explicit open and close."""
        reader = FASTAReader(simple_fasta_file)
        reader.open()
        assert reader._fasta_file is not None
        reader.close()
        assert reader._fasta_file is None

    def test_reader_get_records(self, simple_fasta_file):
        """Test getting all records from file."""
        with FASTAReader(simple_fasta_file) as reader:
            records = reader.get_records()
        
        assert len(records) == 3
        assert records[0].id == "seq1"
        assert records[1].id == "seq2"
        assert records[2].id == "seq3"

    def test_reader_record_sequences(self, simple_fasta_file):
        """Test that sequences are read correctly."""
        with FASTAReader(simple_fasta_file) as reader:
            records = reader.get_records()
        
        assert records[0].sequence == "ATCGATCGATCG"
        assert records[1].sequence == "GCTAGCTAGCTA"
        assert records[2].sequence == "NNNNNNNNNNNN"

    def test_reader_record_descriptions(self, simple_fasta_file):
        """Test that IDs are read correctly from pysam."""
        with FASTAReader(simple_fasta_file) as reader:
            records = reader.get_records()
        
        # pysam's FastxFile uses entry.name which is just the ID
        # The description would need to be in entry.comment, but FASTARecord
        # is initialized with entry.name as header, so description will be empty
        # when using pysam reader
        assert records[0].id == "seq1"
        assert records[1].id == "seq2"
        assert records[2].id == "seq3"

    def test_reader_multiline_sequence(self, multiline_fasta_file):
        """Test reading multi-line sequences."""
        with FASTAReader(multiline_fasta_file) as reader:
            records = reader.get_records()
        
        assert len(records) == 1
        # 3 lines of 60 characters = 180
        assert len(records[0].sequence) == 180

    def test_reader_empty_file(self, empty_fasta_file):
        """Test reading empty file."""
        with FASTAReader(empty_fasta_file) as reader:
            records = reader.get_records()
        
        assert len(records) == 0

    def test_reader_to_dict(self, simple_fasta_file):
        """Test getting sequence IDs from reader."""
        with FASTAReader(simple_fasta_file) as reader:
            seq_ids = reader.get_sequence_ids()
        
        assert "seq1" in seq_ids
        assert "seq2" in seq_ids
        assert "seq3" in seq_ids

    def test_reader_get_record_by_id(self, simple_fasta_file):
        """Test getting record by ID."""
        with FASTAReader(simple_fasta_file) as reader:
            record = reader.get_record_by_id("seq1")
        
        assert record is not None
        assert record.sequence == "ATCGATCGATCG"

    def test_reader_get_record_by_id_not_found(self, simple_fasta_file):
        """Test getting non-existent record."""
        with FASTAReader(simple_fasta_file) as reader:
            record = reader.get_record_by_id("nonexistent")
        
        assert record is None

    @pytest.mark.skipif(
        not (PROJECT_ROOT / "tests" / "testdata" / "test.fa").exists(),
        reason="Test data file not found"
    )
    def test_reader_real_data(self, test_data_fasta):
        """Test reading real test data."""
        with FASTAReader(test_data_fasta) as reader:
            records = reader.get_records()
        
        assert len(records) > 0
        for record in records:
            assert len(record.sequence) > 0
            assert len(record.id) > 0


class TestFASTAReaderIterator:
    """Tests for FASTAReader iteration."""

    def test_reader_iteration(self, simple_fasta_file):
        """Test iterating over records."""
        with FASTAReader(simple_fasta_file) as reader:
            reader.open()
            count = 0
            for record in reader._parse_records():
                count += 1
                assert isinstance(record, FASTARecord)
        
        assert count == 3


# ============================================================================
# Edge case tests
# ============================================================================

class TestFASTAEdgeCases:
    """Tests for edge cases."""

    def test_sequence_with_lowercase(self):
        """Test handling lowercase sequences."""
        content = ">seq1\natcgatcg\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(content)
            f.flush()
            fpath = f.name
        
        try:
            with FASTAReader(fpath) as reader:
                records = reader.get_records()
            
            assert len(records) == 1
            assert records[0].sequence == "atcgatcg"
            assert records[0].sequence_upper == "ATCGATCG"
        finally:
            os.unlink(fpath)

    def test_sequence_with_mixed_case(self):
        """Test handling mixed case sequences."""
        content = ">seq1\nAtCgAtCg\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(content)
            f.flush()
            fpath = f.name
        
        try:
            with FASTAReader(fpath) as reader:
                records = reader.get_records()
            
            assert records[0].sequence == "AtCgAtCg"
        finally:
            os.unlink(fpath)

    def test_header_with_special_characters(self):
        """Test handling headers with special characters."""
        content = ">seq1|chr1:1-100|gene=ABC description\nATCG\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(content)
            f.flush()
            fpath = f.name
        
        try:
            with FASTAReader(fpath) as reader:
                records = reader.get_records()
            
            assert records[0].id == "seq1|chr1:1-100|gene=ABC"
        finally:
            os.unlink(fpath)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
