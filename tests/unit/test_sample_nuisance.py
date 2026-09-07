"""The sample-level endpoint nuisance fit, in the dimension the rest of the
suite does not cover: WHICH reads are admitted to it, and which kernel each of
them is then scored under. The fit and the channel share one admission RULE, so
they cannot disagree about which polyA calls carry 3' evidence -- they do not see
identical read sets, since exclusion and availability intervene."""
from __future__ import annotations

class TestNuisanceAcceptsAdapter:
    """The fit and the channel must share one admission RULE, so ``accept_adapter``
    widens the fit set the same way it widens the channel. Not the same read SET:
    interval exclusion and fit eligibility still intervene."""

    def _inputs(self, qc):
        polya = {f"r{i}": (None if qc == "ADAPTER" else 40.0, qc) for i in range(200)}
        ends = {f"r{i}": ("c", "+", 1000 + (i % 5)) for i in range(200)}
        claimed = {("c", "+"): [1000]}
        return polya, ends, claimed

    def test_adapter_reads_reach_the_distances_the_fit_consumes(self):
        """Not just the counter: the ADAPTER reads must appear in `tes` (the
        annotated-distance path) AND in the per-locus ends the annotation-free
        split-half path builds its modes from."""
        from fin.analysis.sample_nuisance import sample_nuisance_evidence

        polya, ends, claimed = self._inputs("ADAPTER")
        off = sample_nuisance_evidence(polya, ends, claimed_tes=claimed)
        on = sample_nuisance_evidence(polya, ends, claimed_tes=claimed,
                                      accept_adapter=True)
        assert off.n_pass == 0 and off.tes == () and off.ends_by_locus == {}
        assert on.n_pass == 200
        assert {rid for rid, _d in on.tes} == set(polya)
        assert sorted(on.ends_by_locus[("c", "+")]) == sorted(
            e for _c, _s, e in ends.values())
        assert set(on.ids_by_locus[("c", "+")]) == set(polya)
        assert off.n_scorable == on.n_scorable == 200

    def test_pass_population_and_fitted_parameters_are_unchanged_by_the_flag(self):
        from fin.analysis.sample_nuisance import (fit_from_evidence,
                                                  sample_nuisance_evidence)

        polya, ends, claimed = self._inputs("PASS")
        off = sample_nuisance_evidence(polya, ends, claimed_tes=claimed)
        on = sample_nuisance_evidence(polya, ends, claimed_tes=claimed,
                                      accept_adapter=True)
        assert off.n_pass == on.n_pass == 200
        assert off.tes == on.tes and off.ends_by_locus == on.ends_by_locus
        a, b = fit_from_evidence(off).endpoint, fit_from_evidence(on).endpoint
        assert (a.b_narrow, a.b_broad, a.eta, a.n_fit) == (b.b_narrow, b.b_broad,
                                                           b.eta, b.n_fit)


class TestQcConditionalKernels:
    """One endpoint kernel per krill QC class. The point is that HOW FAR a class
    is trusted stops being a constant: it is estimated from the sample's own
    ends, so a class whose ends do not cluster abstains on its own fit rather
    than by appearing on a hand-written list."""

    def _evidence(self, spec, claimed=1000, **kw):
        """spec: qc -> list of 3' ends. Distances are to a single claimed TES."""
        from fin.analysis.sample_nuisance import sample_nuisance_evidence

        polya, ends = {}, {}
        i = 0
        for qc, positions in spec.items():
            for p in positions:
                rid = f"{qc}_{i}"; i += 1
                polya[rid] = (80.0, qc)
                ends[rid] = ("c", "+", p)
        return sample_nuisance_evidence(
            polya, ends, claimed_tes={("c", "+"): [claimed]},
            qc_conditional=True, **kw)

    def test_a_class_whose_ends_do_not_cluster_abstains(self):
        """Abstention is the `available` rule -- the narrow SCALE and the count."""
        from fin.analysis.sample_nuisance import NARROW_SCALE_MAX_BP, fit_from_evidence

        tight = [1000 + (i % 3) - 1 for i in range(300)]
        scattered = [1000 + (i * 37 % 4001) - 2000 for i in range(300)]
        model = fit_from_evidence(self._evidence({"PASS": tight, "NOREGION": scattered}))
        assert model.endpoint_for("PASS") is not None
        assert model.endpoint_by_qc["NOREGION"].b_narrow > NARROW_SCALE_MAX_BP
        assert model.endpoint_for("NOREGION") is None
        assert model.endpoint_for(None) is None
        assert model.endpoint_for("MADE_UP") is None

    def test_too_few_ends_is_also_abstention(self):
        from fin.analysis.sample_nuisance import MIN_KERNEL_FIT, fit_from_evidence

        tight = [1000 + (i % 3) - 1 for i in range(300)]
        few = [1000 + (i % 3) - 1 for i in range(MIN_KERNEL_FIT - 1)]
        model = fit_from_evidence(self._evidence({"PASS": tight, "ADAPTER": few}))
        assert model.endpoint_by_qc["ADAPTER"].n_fit < MIN_KERNEL_FIT
        assert model.endpoint_for("ADAPTER") is None

    def test_the_released_eta_cap_calibrates_a_narrow_but_mostly_wrong_class(self):
        """A class can fit a genuinely narrow component and still be mostly
        read-through (SIRV4 NOREGION: b_narrow 1.8 bp, 34.5% of ends near a TES).
        The pooled cap of 0.5 forces it to claim half its mass is jitter; the
        per-class fit must be free to say otherwise, and that changes the WEIGHT
        a read of that class carries, not whether it is admitted."""
        from fin.analysis.sample_nuisance import (ETA_MAX, ETA_MAX_PER_CLASS,
                                                  fit_endpoint_model)

        # 45% at the mode, the rest spread over +/-200 bp: a sub-bp narrow
        # component that still explains under half the class
        d = [(i % 3) - 1 if i % 20 < 9 else ((i * 53) % 401) - 200 for i in range(600)]
        capped = fit_endpoint_model(d, eta_max=ETA_MAX)
        free = fit_endpoint_model(d, eta_max=ETA_MAX_PER_CLASS)
        assert capped.eta == ETA_MAX and free.eta > ETA_MAX
        assert capped.available and free.available          # both still admitted
        # the released fit trusts a read of this class less: a 100 bp discrepancy
        # costs it less relative to a perfect one
        assert (free.loglik(0) - free.loglik(100)) < (capped.loglik(0)
                                                      - capped.loglik(100))

    def test_an_informative_class_fits_the_same_eta_under_either_cap(self):
        from fin.analysis.sample_nuisance import (ETA_MAX, ETA_MAX_PER_CLASS,
                                                  fit_endpoint_model)

        d = [(i % 3) - 1 if i % 8 else ((i * 29) % 601) - 300 for i in range(600)]
        a = fit_endpoint_model(d, eta_max=ETA_MAX)
        b = fit_endpoint_model(d, eta_max=ETA_MAX_PER_CLASS)
        assert a.eta < ETA_MAX
        assert round(a.eta, 3) == round(b.eta, 3)
        assert round(a.b_narrow, 6) == round(b.b_narrow, 6)

    def test_a_noisier_class_is_trusted_less_not_excluded(self):
        from fin.analysis.sample_nuisance import fit_from_evidence

        tight = [1000 + (i % 3) - 1 for i in range(400)]
        loose = [1000 + (i % 3) - 1 if i % 4 else 1000 + 300 + i for i in range(400)]
        model = fit_from_evidence(self._evidence({"PASS": tight, "READ_FAILED": loose}))
        a, b = model.endpoint_for("PASS"), model.endpoint_for("READ_FAILED")
        assert a is not None and b is not None
        assert b.eta > a.eta            # more of the noisy class is read-through
        # the comparison the channel actually makes is WITHIN one read, between
        # two candidates: a 300 bp discrepancy costs the noisy class less
        assert (b.loglik(0) - b.loglik(300)) < (a.loglik(0) - a.loglik(300))

    def test_exclusion_applies_to_the_per_class_fit(self):
        """Leave-one-interval-out must remove a read from its CLASS's kernel too,
        not only from the pooled one."""
        from fin.analysis.sample_nuisance import fit_from_evidence

        ev = self._evidence({"PASS": [1000 + (i % 3) - 1 for i in range(300)]})
        drop = [rid for rid, _d in ev.tes][:100]
        full = fit_from_evidence(ev).endpoint_by_qc["PASS"]
        held = fit_from_evidence(ev, exclude_reads=drop).endpoint_by_qc["PASS"]
        assert full.n_fit == 300 and held.n_fit == 200

    def test_the_annotation_free_path_also_fits_per_class(self):
        """With no claimed TES the distances come from the split-half read modes;
        the read ids must survive that path or no class can be formed."""
        from fin.analysis.sample_nuisance import fit_from_evidence, sample_nuisance_evidence

        polya, ends = {}, {}
        for qc in ("PASS", "ADAPTER"):
            for i in range(300):
                rid = f"{qc}_{i}"
                polya[rid] = (80.0, qc)
                ends[rid] = ("c", "+", 1000 + (i % 3) - 1)
        ev = sample_nuisance_evidence(polya, ends, qc_conditional=True)
        assert ev.tes == ()                       # no annotation: fallback path
        model = fit_from_evidence(ev)
        assert set(model.endpoint_by_qc) == {"PASS", "ADAPTER"}
        assert model.endpoint_for("PASS") is not None

    def test_the_pooled_model_is_untouched_when_not_conditional(self):
        from fin.analysis.sample_nuisance import (fit_endpoint_model, fit_from_evidence,
                                                  sample_nuisance_evidence)

        polya = {f"r{i}": (80.0, "PASS") for i in range(200)}
        ends = {f"r{i}": ("c", "+", 1000 + (i % 5)) for i in range(200)}
        claimed = {("c", "+"): [1000]}
        ev = sample_nuisance_evidence(polya, ends, claimed_tes=claimed)
        plain = fit_from_evidence(ev)
        assert plain.endpoint_by_qc == {} and plain.qc_conditional is False
        assert plain.endpoint_for("PASS") is plain.endpoint
        assert plain.endpoint_for("ANYTHING") is plain.endpoint
        # and the pooled fit is exactly the mixture fitted on those distances
        ref = fit_endpoint_model([d for _rid, d in ev.tes])
        e = plain.endpoint
        assert (e.b_narrow, e.b_broad, e.eta, e.n_fit) == (ref.b_narrow, ref.b_broad,
                                                           ref.eta, ref.n_fit)

    def test_conditional_evidence_with_no_fittable_class_does_not_fall_back(self):
        """An empty class map in CONDITIONAL mode means no class qualified -- not
        permission to use the pooled kernel. The model and the assignment channel
        must agree about that."""
        from fin.analysis.sample_nuisance import fit_from_evidence

        spec = {f"QC{i}": [1000 + (i % 3) - 1] for i in range(60)}
        model = fit_from_evidence(self._evidence(spec))
        assert model.endpoint_by_qc == {} and model.qc_conditional is True
        assert model.endpoint_for("QC0") is None and model.endpoint_for(None) is None

    def test_the_pooled_kernel_keeps_the_name_admitted_population(self):
        """The pooled kernel's consumer is the discovery contradiction test,
        whose threshold is not invariant to eta, so widening admission must not
        move it: it is fitted on the reads the VERDICT-NAME rule admits, and the
        widened classes reach only the per-class kernels."""
        from fin.analysis.sample_nuisance import fit_endpoint_model, fit_from_evidence

        tight = [1000 + (i % 3) - 1 for i in range(300)]
        scattered = [1000 + (i * 37 % 4001) - 2000 for i in range(300)]
        ev = self._evidence({"PASS": tight, "NOREGION": scattered})
        assert set(ev.pooled.qc_by_read) == set()          # the name-based derivation
        assert {rid for rid, _d in ev.pooled.tes} == {rid for rid in ev.qc_by_read
                                                      if rid.startswith("PASS_")}
        model = fit_from_evidence(ev)
        pooled, cls = model.endpoint, model.endpoint_by_qc["PASS"]
        # the pooled kernel IS the PASS kernel, not a fit over the widened set
        assert pooled.n_fit == 300
        assert (pooled.b_narrow, pooled.b_broad) == (cls.b_narrow, cls.b_broad)
        ref = fit_endpoint_model([d for rid, d in ev.tes if rid.startswith("PASS_")])
        assert (pooled.b_narrow, pooled.eta, pooled.n_fit) == (ref.b_narrow, ref.eta,
                                                               ref.n_fit)

    def test_the_released_cap_reaches_the_class_kernels_through_the_fit(self):
        """The wiring, not just the fitter: fit_from_evidence must supply
        ETA_MAX_PER_CLASS, or a class can never report an eta above the pooled
        cap no matter what its ends look like."""
        from fin.analysis.sample_nuisance import (ETA_MAX, ETA_MAX_PER_CLASS,
                                                  fit_endpoint_model, fit_from_evidence)

        binding = [1000 + ((i % 3) - 1 if i % 20 < 9 else ((i * 53) % 401) - 200)
                   for i in range(600)]
        tight = [1000 + (i % 3) - 1 for i in range(300)]
        model = fit_from_evidence(self._evidence({"PASS": tight, "NOREGION": binding}))
        got = model.endpoint_by_qc["NOREGION"]
        assert got.eta > ETA_MAX
        ref = fit_endpoint_model([x - 1000 for x in binding], eta_max=ETA_MAX_PER_CLASS)
        assert round(got.eta, 9) == round(ref.eta, 9)

    def test_a_fitted_model_stays_usable_as_a_dict_key(self):
        """Frozen means callers may key caches on a fitted model; a dict field
        would silently break that."""
        from fin.analysis.sample_nuisance import fit_from_evidence, fit_sample_nuisance

        tight = [1000 + (i % 3) - 1 for i in range(300)]
        conditional = fit_from_evidence(self._evidence({"PASS": tight}))
        assert {conditional: 1}[conditional] == 1
        assert {fit_sample_nuisance({}, {}): 1}                 # the empty case too

    def test_the_mode_is_recorded_not_inferred_from_the_verdicts(self):
        """A conditional derivation over a sample with no usable verdict is still
        a conditional derivation: inferring the mode from `qc_by_read` would
        report pooled mode and re-open the fallback this design closed."""
        from fin.analysis.sample_nuisance import fit_from_evidence, sample_nuisance_evidence

        ev = sample_nuisance_evidence({}, {}, qc_conditional=True)
        assert ev.qc_by_read == {} and ev.qc_conditional is True
        model = fit_from_evidence(ev)
        assert model.qc_conditional is True
        assert model.endpoint_for("PASS") is None

    def test_the_nested_pooled_fit_survives_reordering_and_exclusion(self):
        """Why the pooled kernel is re-DERIVED under the name-based rule instead
        of filtered out of the widened distance list: admitting more reads can
        change which locus is seen first, the distance list is summed in that
        order, and a filtered list is then the same numbers in a different order,
        so the EM lands a few ULPs away. The read ids here are chosen so the two
        derivations really do order the loci differently."""
        from fin.analysis.sample_nuisance import (fit_endpoint_model, fit_from_evidence,
                                                  sample_nuisance_evidence)

        polya, ends = {}, {}
        claimed = {("cA", "+"): [1000], ("cB", "+"): [900000]}
        for i in range(400):                       # locus A: PASS only, sorts 2nd
            polya[f"m{i}"] = (80.0, "PASS")
            ends[f"m{i}"] = ("cA", "+", 1000 + (i * 7 % 13) - 6)
        for i in range(400):                       # locus B: PASS, sorts 3rd
            polya[f"z{i}"] = (80.0, "PASS")
            ends[f"z{i}"] = ("cB", "+", 900000 + (i * 11 % 23) - 11)
        for i in range(400):                       # locus B: NOREGION, sorts 1st
            polya[f"a{i}"] = (-1.0, "NOREGION")
            ends[f"a{i}"] = ("cB", "+", 900000 + (i * 53 % 401) - 200)
        kw = dict(claimed_tes=claimed)
        plain = sample_nuisance_evidence(polya, ends, **kw)
        cond = sample_nuisance_evidence(polya, ends, qc_conditional=True, **kw)
        assert list(plain.ends_by_locus) == [("cA", "+"), ("cB", "+")]
        assert list(cond.ends_by_locus) == [("cB", "+"), ("cA", "+")]

        def fields(m):
            return (m.b_narrow, m.b_broad, m.eta, m.n_fit)

        drop = sorted(polya)[:300]
        for excl in ((), drop):
            a = fit_from_evidence(plain, exclude_reads=excl).endpoint
            b = fit_from_evidence(cond, exclude_reads=excl).endpoint
            assert fields(a) == fields(b)
        # the rejected implementation -- filter the widened list by read id --
        # really is different, so this fixture would catch a regression to it
        keep = {rid for rid in cond.qc_by_read if not rid.startswith("a")}
        filtered = fit_endpoint_model([d for rid, d in cond.tes if rid in keep])
        assert fields(filtered) != fields(fit_from_evidence(plain).endpoint)

        # and the same through the annotation-free split-half path, where the
        # widened ends would otherwise also move the MODES the fit measures from
        p2 = sample_nuisance_evidence(polya, ends)
        c2 = sample_nuisance_evidence(polya, ends, qc_conditional=True)
        for excl in ((), drop):
            a = fit_from_evidence(p2, exclude_reads=excl).endpoint
            b = fit_from_evidence(c2, exclude_reads=excl).endpoint
            assert fields(a) == fields(b) and a.n_fit > 0
