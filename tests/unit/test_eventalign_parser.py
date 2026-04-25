"""Tests for eventalign TSV parser and distance matrix construction."""

import math
import sys
import tempfile
from pathlib import Path


PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fin.scoring.eventalign_parser import (
    ReadCandidateScore,
    build_distance_matrix,
    parse_eventalign_tsv,
)

HEADER = "contig\tposition\treference_kmer\tread_name\tstrand\tevent_index\tevent_level_mean\tevent_stdv\tevent_length\tmodel_kmer\tmodel_mean\tmodel_stdv\tstandardized_level\tstart_idx\tend_idx\tsamples"


def _make_row(contig, position, read_name, event_level_mean, model_mean, model_stdv, start_idx, end_idx):
    return f"{contig}\t{position}\tACGTT\t{read_name}\t+\t0\t{event_level_mean}\t1.0\t0.01\tACGTT\t{model_mean}\t{model_stdv}\t0.0\t{start_idx}\t{end_idx}\tsamples"


class TestParseEventalignTsv:
    """Tests for TSV parsing."""

    def test_basic_parsing(self):
        """Parse a minimal eventalign TSV."""
        lines = [
            HEADER,
            _make_row("tx1", 0, "read1", 100.0, 100.0, 1.0, 0, 100),
            _make_row("tx1", 1, "read1", 101.0, 100.0, 1.0, 100, 200),
            _make_row("tx2", 0, "read1", 90.0, 100.0, 1.0, 0, 100),
            _make_row("tx1", 0, "read2", 99.0, 100.0, 1.0, 50, 150),
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("\n".join(lines) + "\n")
            tsv_path = f.name

        scores = parse_eventalign_tsv(tsv_path, candidate_lengths={"tx1": 100, "tx2": 50})

        # Should have 3 pairs: (read1, tx1), (read1, tx2), (read2, tx1)
        assert len(scores) == 3

        r1_tx1 = [s for s in scores if s.read_name == "read1" and s.candidate_id == "tx1"][0]
        assert r1_tx1.num_events == 2
        assert r1_tx1.coverage == 2 / 100  # 2 unique positions / 100 length
        assert r1_tx1.signal_start_idx == 0
        assert r1_tx1.signal_end_idx == 200

    def test_no_header(self):
        """Parse TSV without header row."""
        lines = [
            _make_row("tx1", 0, "read1", 100.0, 100.0, 1.0, 0, 100),
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("\n".join(lines) + "\n")
            tsv_path = f.name

        scores = parse_eventalign_tsv(tsv_path)
        assert len(scores) == 1

    def test_log_likelihood_perfect_match(self):
        """When event == model, log-likelihood should be near maximum."""
        lines = [
            HEADER,
            _make_row("tx1", 0, "read1", 100.0, 100.0, 1.0, 0, 100),
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("\n".join(lines) + "\n")
            tsv_path = f.name

        scores = parse_eventalign_tsv(tsv_path)
        s = scores[0]
        # z=0, ll = -0.5 * (0 + log(2pi) + 0) = -0.5 * log(2pi) ≈ -0.919
        expected_ll = -0.5 * math.log(2 * math.pi)
        assert abs(s.total_log_likelihood - expected_ll) < 0.01

    def test_rmse(self):
        """RMSE should reflect deviation from model mean."""
        lines = [
            HEADER,
            _make_row("tx1", 0, "read1", 102.0, 100.0, 1.0, 0, 100),
            _make_row("tx1", 1, "read1", 98.0, 100.0, 1.0, 100, 200),
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("\n".join(lines) + "\n")
            tsv_path = f.name

        scores = parse_eventalign_tsv(tsv_path)
        s = scores[0]
        # errors: (2)^2=4, (-2)^2=4, mean=4, sqrt=2
        assert abs(s.event_rmse - 2.0) < 0.01


class TestBuildDistanceMatrix:
    """Tests for distance matrix construction."""

    def test_basic_matrix(self):
        """Build distance matrix from scores."""
        scores = [
            ReadCandidateScore(read_name="r1", candidate_id="tx1", total_log_likelihood=-10.0),
            ReadCandidateScore(read_name="r1", candidate_id="tx2", total_log_likelihood=-20.0),
            ReadCandidateScore(read_name="r2", candidate_id="tx1", total_log_likelihood=-15.0),
        ]

        read_ids = ["r1", "r2"]
        candidate_ids = ["tx1", "tx2"]

        dist = build_distance_matrix(scores, read_ids, candidate_ids)

        assert dist.shape == (2, 2)
        assert dist[0, 0] == 10.0  # -(-10)
        assert dist[0, 1] == 20.0  # -(-20)
        assert dist[1, 0] == 15.0  # -(-15)
        assert dist[1, 1] == 1e6  # missing pair

    def test_missing_reads_get_large_distance(self):
        """Reads/candidates with no scores should have large distance."""
        dist = build_distance_matrix([], ["r1"], ["tx1"])
        assert dist[0, 0] == 1e6
