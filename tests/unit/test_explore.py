"""Unit tests for the explore step (family-local intron-graph path enumeration).

explore_family_paths ADDS assembled full-length candidates that no single read produced, by tiling
read-phased junction edges across reads within one family. It must never drop/rewrite; it only emits
NEW (non-variant) maximal source->sink paths.
"""
from fin.candidates.chain_cluster import ChainFamily
from fin.candidates.explore import ExploredPath, explore_family_paths

# convenient introns, genomically ordered with positive exons between
I1 = (1000, 2000)
I2 = (3000, 4000)
I3 = (5000, 6000)
I4 = (7000, 8000)


def _fam(variant_reads):
    """variant_reads: {chain_tuple: set(read_ids)} -> a ChainFamily."""
    variants = list(variant_reads)
    pool = set()
    for r in variant_reads.values():
        pool |= r
    return ChainFamily(variants=variants, read_pool=pool, variant_reads=dict(variant_reads))


def _chains(paths):
    return sorted(p.chain for p in paths)


class TestCrossReadAssembly:
    def test_truncated_reads_assemble_full_chain(self):
        # read A = [I1,I2], read B = [I2,I3]; NO read has [I1,I2,I3] -> assemble it.
        fam = _fam({(I1, I2): {"a1", "a2"}, (I2, I3): {"b1"}})
        out = explore_family_paths(fam)
        assert _chains(out) == [(I1, I2, I3)]
        p = out[0]
        assert p.min_edge_support == 1          # bottleneck = edge I2->I3 (1 read)
        assert p.n_edges == 2
        assert p.reads == {"a1", "a2", "b1"}     # union of tiling variants' reads

    def test_three_way_tiling(self):
        fam = _fam({(I1, I2): {"a"}, (I2, I3): {"b"}, (I3, I4): {"c"}})
        out = explore_family_paths(fam)
        assert (I1, I2, I3, I4) in _chains(out)   # full 4-intron chain assembled

    def test_single_variant_emits_nothing(self):
        fam = _fam({(I1, I2): {"a"}})
        assert explore_family_paths(fam) == []

    def test_full_read_chain_not_re_emitted(self):
        # one read already spans the whole chain -> it's a variant -> explore adds nothing new.
        fam = _fam({(I1, I2, I3): {"a"}})
        assert explore_family_paths(fam) == []


class TestFrankensteinRecombination:
    # A-B-X-D and A-C-X-E share node X; read-phased edges also permit A-B-X-E / A-C-X-D
    # (each edge observed, whole chain not). Codex flagged this global recombination; we EMIT
    # them (recall-first) -- this test documents/pins that behaviour.
    A = (1000, 2000)
    B = (2500, 3000)
    C = (2500, 3100)
    X = (4000, 5000)
    D = (6000, 7000)
    E = (6000, 7100)

    def test_recombinants_are_emitted(self):
        v1 = (self.A, self.B, self.X, self.D)
        v2 = (self.A, self.C, self.X, self.E)
        fam = _fam({v1: {"r1"}, v2: {"r2"}})
        out = _chains(explore_family_paths(fam))
        recomb1 = (self.A, self.B, self.X, self.E)
        recomb2 = (self.A, self.C, self.X, self.D)
        assert recomb1 in out and recomb2 in out   # both recombinants produced
        assert v1 not in out and v2 not in out      # variants are NOT re-emitted (kept by caller)


class TestBeamAndDeterminism:
    def test_deterministic_across_variant_order(self):
        d1 = {(I1, I2): {"a"}, (I2, I3): {"b"}}
        d2 = {(I2, I3): {"b"}, (I1, I2): {"a"}}
        assert _chains(explore_family_paths(_fam(d1))) == _chains(explore_family_paths(_fam(d2)))

    def test_max_new_caps_output(self):
        A = (1000, 2000)
        B, C = (2500, 3000), (2500, 3100)
        X = (4000, 5000)
        D, E = (6000, 7000), (6000, 7100)
        fam = _fam({(A, B, X, D): {"r1"}, (A, C, X, E): {"r2"}})
        out = explore_family_paths(fam, max_new=1)
        assert len(out) == 1                        # capped

    def test_returns_exploredpath_objects(self):
        fam = _fam({(I1, I2): {"a"}, (I2, I3): {"b"}})
        out = explore_family_paths(fam)
        assert all(isinstance(p, ExploredPath) for p in out)

    def test_variant_sinks_do_not_starve_budget(self):
        # many already-observed full variants PLUS one cross-read tiled chain. The exact variants
        # get discarded downstream, so they must not starve the real assembled chain out of the
        # (small) output budget.
        variants = {
            ((100, 200), (300, 400), (500, 600)): {"v1"},
            ((1000, 1100), (1200, 1300), (1400, 1500)): {"v2"},
            ((2000, 2100), (2200, 2300), (2400, 2500)): {"v3"},
            ((3000, 3100), (3200, 3300), (3400, 3500)): {"v4"},
            ((4000, 4100), (4200, 4300), (4400, 4500)): {"v5"},
        }
        # a family whose reads tile a 3-intron chain that no single read spans
        variants[(I1, I2)] = {"t1"}
        variants[(I2, I3)] = {"t2"}
        out = explore_family_paths(_fam(variants), max_new=1)
        assert (I1, I2, I3) in _chains(out)   # the assembled chain is not starved out

    def test_strongest_assembled_path_wins_small_budget(self):
        # Codex gpt-5.6-sol HIGH repro: several shallow support-1 tiled sinks vs one longer
        # support-10 tiled sink. With max_new=1 the loop must NOT stop early on a completion
        # count -- ranking is bottleneck-first, so the support-10 chain must be the one returned.
        vr = {}
        for i in range(5):                       # five shallow, support-1 assembled sinks
            base = i * 1000
            s, mid, t = (base, base + 10), (base + 20, base + 30), (base + 40, base + 50)
            vr[(s, mid)] = {f"a{i}"}
            vr[(mid, t)] = {f"b{i}"}
        L = ((10000, 10010), (10020, 10030), (10040, 10050), (10060, 10070))
        for j, (a, b) in enumerate(zip(L, L[1:])):   # one longer, support-10 assembled sink
            vr[(a, b)] = {f"L{j}_{k}" for k in range(10)}
        out = explore_family_paths(_fam(vr), max_new=1)
        assert len(out) == 1
        assert out[0].chain == L                 # strongest (support-10) path, not a shallow one
        assert out[0].min_edge_support == 10

    def test_max_new_zero_returns_empty(self):
        fam = _fam({(I1, I2): {"a"}, (I2, I3): {"b"}})
        assert explore_family_paths(fam, max_new=0) == []

    def test_beam1_does_not_drop_high_support_source(self):
        # Codex gpt-5.6-sol MED: the initial source frontier must NOT be beam-truncated -- one-node
        # paths all score (0,0), so a beam=1 cut on sources would deterministically drop a later
        # high-support component. Two disjoint tiled chains; the support-10 one must survive beam=1.
        weak_a, weak_m, weak_t = (100, 200), (300, 400), (500, 600)
        strong = ((9000, 9010), (9020, 9030), (9040, 9050))
        vr = {
            (weak_a, weak_m): {"wa"}, (weak_m, weak_t): {"wb"},          # support-1 tiled chain
            (strong[0], strong[1]): {f"s{k}" for k in range(10)},        # support-10 ...
            (strong[1], strong[2]): {f"t{k}" for k in range(10)},        # ... tiled chain
        }
        out = explore_family_paths(_fam(vr), beam_width=1, max_new=1)
        assert out and out[0].chain == strong    # strongest survives even at beam_width=1


class TestGeometryGuard:
    def test_no_edge_across_overlapping_introns(self):
        # a crafted (invalid) variant whose consecutive introns overlap must not form an edge,
        # so no assembled path bridges them.
        bad = ((1000, 4000), (3000, 5000))          # second donor 3000 <= first acceptor 4000
        other = ((3000, 5000), (6000, 7000))
        fam = _fam({bad: {"a"}, other: {"b"}})
        out = _chains(explore_family_paths(fam))
        # (1000,4000) has no valid outgoing edge -> cannot bridge into (6000,7000)
        assert ((1000, 4000), (3000, 5000), (6000, 7000)) not in out
