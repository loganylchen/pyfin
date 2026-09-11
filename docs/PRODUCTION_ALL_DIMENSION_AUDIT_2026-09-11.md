# Where is this tool from optimal, dimension by dimension?

An audit rather than a lever. The point is to stop choosing what to work on
from memory and instead measure every axis the tool is judged on, against the
best available alternative on that axis, and say plainly which are settled and
which are open. Everything below is measured on runs on disk, under the
expression contract, unless marked otherwise.

## Summary

| # | dimension | status | best alternative | gap |
| --- | --- | --- | --- | --- |
| 1 | structure recall, annotated | **ahead or level** | stringtie3 | +1.4 mean |
| 2 | structure precision, annotated | **ahead** | stringtie3 | +5.8 mean |
| 3 | robustness to corrupt annotation | **far ahead** | stringtie3 | +10.5 (jitter) |
| 4 | robustness to sparse annotation | **ahead** | stringtie3 | +3.6 (p50) |
| 5 | **de novo (p00)** | **BEHIND** | stringtie3 | **−0.30** |
| 6 | transcript ends | not a gap in the matched set | — | ends exact |
| 7 | quantification | bottlenecked by detection | NanoCount oracle | — |
| 8 | signal usage | genuinely discriminative | nothing published | — |

One structural dimension is open: **de novo**. Everything else is level or
ahead, and dimension 7 is downstream of 1.

## 1–4. Structure, against the one real peer

Real human (H9 rep2 run2, GENCODE v44, chr1/19/21, 7 conditions), after the
attestation exemption landed:

| cond | dSn | dPr | dF1 |
| --- | ---: | ---: | ---: |
| `c_jitter50_10bp` | +6.9 | +20.8 | **+10.48** |
| `p50` | +2.8 | +4.3 | +3.56 |
| `c_merge50` | +0.8 | +3.4 | +1.44 |
| `c_trunc5p50` | −0.2 | +5.5 | +1.09 |
| `full` | +0.1 | +3.3 | +0.83 |
| `c_trunc3p50` | −0.2 | +3.7 | +0.65 |
| `p00` | −0.2 | −0.6 | **−0.30** |
| **mean** | **+1.4** | **+5.8** | **+2.54** |

Ahead on 6 of 7, and the shape has changed: before attestation we were behind
on 5 of 7 with a 1.8–3.0 pp recall deficit everywhere. Recall is now level
(−0.2 to +0.1) on the intact and truncated cells while precision is +3.3 to
+5.5, and the lead widens exactly where the annotation is degraded.

## 5. De novo is the one open structural dimension

`p00` (no annotation at all) is −0.2 Sn / −0.6 Pr / −0.30 F1. It is also the
one cell where the attestation exemption is **exactly inert** — an
annotation-conditioned exemption has nothing to readmit when there is no
annotation — so the lever that fixed the other six cannot touch this one.

## 6. Transcript ends are not a gap in the matched set

For models that match a truth transcript, both arms place both ends exactly:

| arm | n | 5' median offset | ≤25 bp | 3' median | ≤25 bp |
| --- | ---: | ---: | ---: | ---: | ---: |
| before attestation | 1903 | 0 | 100% | 0 | 100% |
| after | 2105 | 0 | 100% | 0 | 100% |

The `.ext` read-extent population — models whose ends were replaced by a read
extent — is a separate set that does not receive `=` matches, and was measured
earlier: 86% of its 5' claims do not reproduce in an independent library. The
lever for it exists and stays default-off. Attestation is neutral on this axis.

A methodological note worth keeping: `gffcompare --strict-match` returns **zero**
matches on real dRNA at any `-e`, because it demands exact terminal coordinates
and 5' degradation guarantees they differ. That is the wrong instrument for
this question, not a 0% score; the offset distribution above is the right one.

## 7. Quantification is downstream of detection

Measured earlier on Sequin against known input molarity: the apparent gap to
stringtie3 was a reporting convention (a length-divided `tpm` for a protocol
that is not length-proportional), and the residual over-estimate on long
transcripts is absorbed abundance from **undetected same-gene transcripts**
(r = 0.341; −0.167 clean vs +1.832 with missed siblings; the NanoCount oracle
does not show it). So this axis improves when detection improves — which
dimension 1 just did — and has no independent lever worth pulling first.

The literature agrees on the convention point independently: oarfish,
TranSigner, Bambu and Isosceles all decline length normalisation for
end-anchored long-read protocols, because with no fragmentation
E[reads] ∝ θ, not θ × L.

## 8. Signal is doing real work, and nothing published does this

Measured this session with a signal-free control (same tie population, same
abstentions, same thresholded margin, only the winning structure randomised):

| grid | tie-break worth | of which signal | concentration alone |
| --- | ---: | ---: | --- |
| SIRV4 | −0.93 | −0.75 (81%) | +0.18, CI straddles 0 |
| HEYA8 | −0.76 | −0.59 (78%) | +0.17, CI straddles 0 |
| human | −0.27 | −0.14 (67% of Pr) | +0.13 |

A literature sweep found **no method operating in squiggle space for isoform
assignment at all** — NanoCount, Bambu, LIQA, oarfish, TranSigner and the EM
family all work post-alignment and resolve ambiguity by fractional assignment;
signal-level direct-RNA work is concentrated on modification calling. So this
is a real if small differentiator, worth −0.14 F1 on human.

Per-candidate-base normalisation, the one piece of the NanoSplicer prescription
production does not implement, was tested and is a **null** (SIRV4 +0.08 with a
CI including zero; HEYA8 −0.02 with the family means flipping sign).

## A premise from the literature, tested and mostly not confirmed

minimap2's `--splice-flank=yes` (the default, and what every grid here was
aligned with) scores a GT donor higher when the next base is A or G, a prior
taken from human. Its own cookbook says SIRV controls follow that trend only
~60% of the time and recommends `--splice-flank=no` when benchmarking on SIRVs.
If true in consequence, part of the junction wobble this project has spent
years engineering against would be an aligner artefact specific to the grids
the machinery was tuned on.

The premise holds on our references:

| grid | GT donors followed by A/G |
| --- | ---: |
| human | **93.9%** (prior holds, better than the documented 91%) |
| SIRV4 | **62.3%** (violated) |
| HEYA8 | **50.9%** (essentially random) |

The consequence largely does not. Splitting read junctions by whether the
junction they land on satisfies the prior, and asking how often they are placed
exactly:

| grid | A/G exact | C/T exact | gap |
| --- | ---: | ---: | ---: |
| human | 98.6% | 95.2% | **+3.4 pp** |
| SIRV4 | 96.0% | 93.0% | **+3.0 pp** |
| HEYA8 | 95.1% | 94.8% | **+0.3 pp** |

The per-junction penalty is about the same everywhere, and is *smallest* on the
grid that violates the prior most. So the story "the grids that break the
assumption pay for it" is wrong; what differs between grids is the *proportion*
of junctions in the violating class (6% / 38% / 49%), not the damage each one
takes. Multiplying through, the wobble attributable to this is ≈0.2 pp on
human, ≈1.1 pp on SIRV4 and ≈0.2 pp on HEYA8 — real on SIRV4, roughly five
times human, and still about one percentage point of junction placement.

Not worth a realignment campaign. Worth knowing, and worth the `--splice-flank=no`
flag on any *future* SIRV alignment. A separate check passed cleanly: no grid
uses `--junc-bed`, so there is no annotation leakage through the aligner under
the annotation-degradation contract.

## Where the remaining recall actually dies, re-measured

The audit that chose the last lever said `m2_support_gate` owned 62% of the
expressed truth. Re-run against the **new default**, on the same cell:

| owner of the loss | n | share |
| --- | ---: | ---: |
| `m2_support_gate` | 3945 | **81.8%** |
| survived selection, lost later | 438 | 9.1% |
| `cluster_recheck` | 346 | 7.2% |
| never a candidate (generation) | 53 | 1.1% |
| `extent_alternative` | 41 | 0.9% |

Generation is **not** the bottleneck on annotated human data — 1.1%. The same
gate still owns 81.8% of what is left, because K=3 only reaches the tier with
three or more compatible reads. The next tier down (suffix 1–2 reads, 2160
candidates) is 59% truth against an F1 break-even of 21.4%.

So the next question is the **operating point**, not a new mechanism: an arm
each at K=2 and K=1. On SIRV4 those are **byte-identical to K=3** — the depth
threshold never binds on a high-depth spike-in, only the veto does — which
means lowering K carries zero spike-in risk and the human arms decide it alone.
Those are in flight.
