"""Unit tests for the de-novo intron graph (consensus clustering + read-adjacency
graph + unambiguous maximal-chain extension)."""
from __future__ import annotations

from fin.candidates.junction_graph import (
    assemble_chains,
    build_graph,
    cluster_junctions,
)

J1, J2, J3 = (100, 200), (300, 400), (500, 600)


class TestClustering:
    def test_wobbles_map_to_highest_count_consensus(self):
        # J2 seen 3x, a ±2bp wobble seen 1x -> wobble maps to J2 (the consensus).
        chains = [(J2,), (J2,), (J2,), ((302, 398),)]
        cons = cluster_junctions(chains, tol=6)
        assert cons[J2] == J2
        assert cons[(302, 398)] == J2          # absorbed into the mode

    def test_beyond_tol_separate(self):
        chains = [(J2,), (J2,), ((312, 388),)]
        cons = cluster_junctions(chains, tol=6)
        assert cons[(312, 388)] == (312, 388)  # 12bp -> its own centroid


class TestGraph:
    def test_edge_counts(self):
        cc = [(J1, J2), (J1, J2, J3), (J2, J3)]
        nodes, edges = build_graph(cc)
        assert edges[(J1, J2)] == 2   # chains 0 and 1
        assert edges[(J2, J3)] == 2   # chains 1 and 2
        assert nodes[J2] == 3


class TestAssembly:
    def test_truncated_reads_assemble_to_full_chain(self):
        # Read A = 5' half (J1,J2), read B = 3' half (J2,J3), read C = full.
        # Every read should extend to the full chain (J1,J2,J3).
        chains = [(J1, J2), (J2, J3), (J1, J2, J3)]
        _, ext = assemble_chains(chains, tol=6, min_edge_reads=1)
        assert ext[0] == (J1, J2, J3)   # A extended right through J2->J3
        assert ext[1] == (J1, J2, J3)   # B extended left through J1->J2
        assert ext[2] == (J1, J2, J3)

    def test_branch_point_stops_extension(self):
        # J2 has TWO supported successors (J3 and J3b) -> ambiguous -> a read
        # ending at J2 stays truncated (no fabricated combo).
        J3b = (700, 800)
        chains = [(J1, J2), (J2, J3), (J2, J3b)]
        _, ext = assemble_chains(chains, tol=6, min_edge_reads=1)
        assert ext[0] == (J1, J2)       # cannot pick J3 vs J3b -> stops at J2
        assert ext[1] == (J1, J2, J3)   # B: left through J1->J2 (J2 unique pred)
        assert ext[2] == (J1, J2, J3b)

    def test_unsupported_edge_not_followed_for_extension(self):
        # J1->J2 has 2 reads (supported at min_edge_reads=2); J2->J3 has 1 (pruned).
        chains = [(J1, J2), (J1, J2), (J2, J3)]
        _, ext = assemble_chains(chains, tol=6, min_edge_reads=2)
        # read (J1,J2) must NOT extend right through the pruned J2->J3 edge:
        assert ext[0] == (J1, J2)
        # read (J2,J3) keeps its own junctions but extends LEFT through the
        # supported J1->J2 (a 3'-partial getting its 5' completed):
        assert ext[2] == (J1, J2, J3)

    def test_wobbled_partials_assemble(self):
        # Same as full-assembly but read B's J2 is a 2bp wobble -> consensus
        # clustering unifies it so the graph still connects.
        chains = [(J1, J2), ((302, 398), J3), (J1, J2, J3)]
        _, ext = assemble_chains(chains, tol=6, min_edge_reads=1)
        assert ext[0] == (J1, J2, J3)
        assert ext[1] == (J1, J2, J3)   # wobbled J2 clustered to J2, then extended


class TestTssBrake:
    # A genuine short isoform (J2,J3) whose reads all START (5'-end) at one TSS,
    # contained inside a longer transcript (J1,J2,J3). Without the brake the short
    # reads extend left through the unique J1->J2 edge and the short isoform is
    # merged away; with the brake, J2 is a real TSS so they stay (J2,J3).
    def test_plus_strand_tss_peak_blocks_5p_extension(self):
        chains = [(J2, J3), (J2, J3), (J2, J3), (J2, J3), (J1, J2, J3)]
        fp = [250, 250, 252, 248, 150]   # 4 short reads pile at ~250 (TSS peak)
        _, ext = assemble_chains(chains, tol=6, min_edge_reads=1,
                                 five_prime_pos=fp, strand="+")
        assert ext[0] == (J2, J3)        # brake kept the short isoform
        assert ext[4] == (J1, J2, J3)    # long isoform still maximal

    def test_brake_off_merges_the_short_isoform(self):
        # Same data, no five_prime_pos -> byte-identical un-braked behaviour.
        chains = [(J2, J3), (J2, J3), (J2, J3), (J2, J3), (J1, J2, J3)]
        _, ext = assemble_chains(chains, tol=6, min_edge_reads=1)
        assert ext[0] == (J1, J2, J3)    # short reads extended left -> merged

    def test_scattered_5p_ends_do_not_brake(self):
        # Reads starting at J2 but with SCATTERED 5'-ends (degradation, no peak)
        # must NOT brake -> they still extend (this is a true truncation).
        chains = [(J2, J3), (J2, J3), (J2, J3), (J2, J3), (J1, J2, J3)]
        fp = [250, 190, 275, 220, 150]   # spread > tol, no >=40% pile
        _, ext = assemble_chains(chains, tol=6, min_edge_reads=1,
                                 five_prime_pos=fp, strand="+")
        assert ext[0] == (J1, J2, J3)    # no TSS peak -> extends as before

    def test_minus_strand_tss_peak_blocks_5p_extension(self):
        # On '-' strand the 5' end is the LARGE coord (reference_end) and 5'-ward
        # extension is to the RIGHT. Short isoform (J1,J2) with reads piling their
        # 5'-end just past J2; brake blocks the J2->J3 (rightward) extension.
        chains = [(J1, J2), (J1, J2), (J1, J2), (J1, J2), (J1, J2, J3)]
        fp = [450, 450, 452, 448, 650]   # short reads' 5'-end (ref_end) pile ~450
        _, ext = assemble_chains(chains, tol=6, min_edge_reads=1,
                                 five_prime_pos=fp, strand="-")
        assert ext[0] == (J1, J2)        # brake kept the short isoform
        assert ext[4] == (J1, J2, J3)
