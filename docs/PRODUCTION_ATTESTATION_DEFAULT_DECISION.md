# Should the support-gate attestation exemption become a default?

Pre-registered before any arm finished. The promotion rule and the kill test
below are fixed in advance precisely because I already believe this lever
works, which is the condition under which a decision rule is worth writing
down.

## The target, quantified

Real human (SGNex H9 rep2 run2, GENCODE v44, chr1/19/21, expression contract:
truth = 7398 annotated transcripts with NanoCount TPM>0, gffcompare `=`)
against stringtie3, the one real peer:

| cond | pyfin Sn/Pr/F1 | stringtie3 | dF1 | dSn | dPr | Sn needed to tie |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| `full` | 31.2/67.9/42.75 | 33.9/65.5/44.68 | −1.92 | −2.7 | +2.4 | **+2.1** |
| `c_trunc5p50` | 30.5/66.8/41.88 | 33.5/62.3/43.57 | −1.69 | −3.0 | +4.5 | +1.8 |
| `c_trunc3p50` | 30.9/66.5/42.19 | 32.7/63.1/43.08 | −0.88 | −1.8 | +3.4 | +1.0 |
| `p00` | 15.8/56.2/24.67 | 16.0/56.8/24.97 | −0.30 | −0.2 | −0.6 | +0.2 |
| `c_merge50` | 27.1/61.3/37.58 | 28.0/57.6/37.68 | −0.10 | −0.9 | +3.7 | +0.1 |
| `p50` | 24.4/62.6/35.11 | 22.9/59.1/33.01 | +2.10 | +1.5 | +3.5 | — |
| `c_jitter50_10bp` | 16.3/41.9/23.47 | 9.5/20.9/13.06 | +10.41 | +6.8 | +21.0 | — |

Behind on 5 of 7. On four of those five the shape is identical — they recall
1.8–3.0 pp more, we are 2.4–4.5 pp more precise — so **+2.1 pp of recall at
constant precision closes every losing cell**. (Codex's qualification, which I
had glossed: `p00` also loses precision, so it is not purely a recall story,
and `c_merge50` needs only +0.1.)

## The correction that reframes the trade

I had been looking for a recall lever that costs no precision. That is the
wrong bar. Under ordinary additive TP/FP accounting, adding *n* predictions
that recover *t* previously-unmatched truths improves F1 whenever

    t / n  >  F1_old / 2

At `full` that break-even marginal precision is **≈21.4%**, not 67.9%. So a
candidate pool that is 77.5% truth is not marginally acceptable — it is far
above the bar, and the question is only how much of its purity survives as
*marginal* precision (models matching truths we already recover do not count,
and retention changes downstream competition).

## Why this lever and not another

The selection audit attributes **475 of the 554** transcripts stringtie3
matched and we did not on `full`, and **4566 of all 7398** expressed truth
transcripts, to `m2_support_gate`. None of the 554 dies in discovery. A later
abundance floor cannot resurrect a candidate the gate already dropped, so the
gate is the causal bottleneck for the annotated gap. Gate attribution is an
**opportunity ceiling, not a predicted recovery count** — exemption still has
to turn into distinct final truth matches.

Tiering the gate-dropped annotated candidates by read evidence ("suffix" = the
read's intron chain is a strand-aware contiguous 3' suffix of the candidate's,
which is exactly what dRNA 5'-degradation produces):

| tier | n | % truth |
| --- | ---: | ---: |
| exact reads ≥ 3 | 122 | 75% |
| exact 2 / 1 | 71 / 347 | 65% / 45% |
| no exact, **suffix ≥ 3** | **1367** | **77.5%** |
| no exact, suffix 1–2 | 2160 | 59% |
| no chain read at all | 15929 | 12% |

The suffix-only tier is 11× the exact tier at the same purity. Codex confirmed
from the code that at K=3 the exemption **already admits it**: three qualifying
suffix reads suffice with zero exact-chain reads, and there is no hidden
requirement for a complete-chain read
(`fin/candidates/discovery.py:203,225,235,237,246`). So the lever is aimed at
the right pool and the first move is to measure the funnel, not to widen the
rule.

Two corrections to my own description, from Codex reading the code:

* The unjittered veto compares candidate structure against **observed read
  junctions**, not against another candidate's junctions
  (`fin/scoring/m2_junction_nll.py:521,532`).
* Attestation only removes a candidate from the **gate's** drop set. It still
  faces the separate GTF abundance floor
  (`fin/pipeline/selection.py:393,847`). If the arms show candidates merely
  moving from `m2_support_gate` to `min_abundance`, the next experiment is
  attestation × GTF floor, not a novel-floor sweep.

## The failure hypothesis to falsify

Not "precision might drop". The specific mechanism:

> Attestation improves an annotation-conditioned expression score by retaining
> multiple plausible annotated *extensions* of the same shared 3' suffix —
> crediting every compatible isoform rather than recovering independently
> supported ones — and therefore becomes harmful where the annotation is
> incomplete or corrupted.

The same three suffix reads can attest several nested annotated models whose
distinguishing upstream structure was never observed. That is deliberate in
the implementation and is covered by its tests
(`fin/candidates/discovery.py:181`, `tests/unit/test_support_gate_attest.py:69`).
Raising K does **not** cure it: a deeply covered shared suffix credits every
compatible extension no matter how deep.

## Promotion rule, fixed before the numbers

K = 3, veto = `novel`, frozen. Compared against **production**, not against
peer win-counts.

1. Positive mean ΔF1 **separately** on human, HEYA8 and SIRV4. A losing grid
   is not pooled away.
2. **Reject** if any pre-specified corruption family has mean ΔF1 ≤ **−1.0 pp**
   within a grid. Smaller negatives are reported, not ignored; the 1 pp
   tolerance is a policy choice, not a significance test.
3. Positive mean ΔF1 across the three human target cells `full`,
   `c_trunc5p50`, `c_trunc3p50` — otherwise it has not solved the stated
   problem, whatever it did elsewhere.
4. Inspect the suffix-only additions by how many models share their supporting
   reads and whether the distinguishing upstream junctions are observed. If the
   high-sharing additions are net F1-negative, the blanket exemption has failed
   in exactly its vulnerable regime.

156 and 130 cells are not independent biological replicates; effects are
summarised paired, by library and by condition family. Passing justifies
promotion **within the tested operating domain** only — it does not establish
that every exempted isoform is independently expressed.

## Arms

One clean snapshot (HEAD `735cbf8`, `_snap/att_dec`), one lever
(`--m2-support-gate-attest-reads 3`), default control re-run on the same
snapshot for both spike-in grids.

| grid | default | attested |
| --- | --- | --- |
| SIRV4 (156) | `adef` 211946 | `aatt` 211947 |
| HEYA8 (130) | `adef` 211948 | `aatt` 211949 |
| human (7) | `m2base` (verified byte-identical at default) | `aatt` 211950 |


---

# Outcome: promoted (commit `17fec28`)

All four rules passed. The lever is now `real-drna`'s default at K=3.

| grid | dSn | dPr | dF1 | 95% CI |
| --- | ---: | ---: | ---: | --- |
| human GENCODE (7 conds) | **+1.47** | **+0.36** | **+1.45** | one library, no honest CI |
| SIRV4 (156 / 6 libs) | +0.18 | −0.15 | +0.12 | [−0.05, +0.32] |
| HEYA8 (130 / 5) | +0.04 | −0.03 | +0.04 | [+0.01, +0.09] |

| cond | before | after | vs stringtie3 |
| --- | --- | --- | ---: |
| `full` | 31.2/67.9/42.75 | **34.0/68.8/45.51** | +0.83 |
| `c_trunc5p50` | 30.5/66.8/41.88 | **33.3/67.8/44.66** | +1.09 |
| `c_trunc3p50` | 30.9/66.5/42.19 | **32.5/66.8/43.73** | +0.65 |
| `c_merge50` | 27.1/61.3/37.58 | 28.8/61.0/39.13 | +1.44 |
| `p50` | 24.4/62.6/35.11 | 25.7/63.4/36.57 | +3.56 |
| `c_jitter50_10bp` | 16.3/41.9/23.47 | 16.4/41.7/23.54 | +10.48 |
| `p00` | 15.8/56.2/24.67 | *inert* | −0.30 |

**Cells beating stringtie3: 2/7 → 6/7.** Rule 3 (the three target cells)
+2.36. Rule 2: no corruption family negative in the mean on either spike-in
grid; worst single cells −0.43 (SIRV4) and −0.16 (HEYA8) against the −1.0
rejection line. `p00` is exactly inert, which is the right behaviour — with no
annotation an annotation-conditioned exemption has nothing to readmit, and it
is also the one cell still behind.

**Precision rose.** 194 distinct newly-matched truth transcripts from 232 new
models = 84% marginal precision against the 21.4% break-even. The lever did not
spend the precision surplus; it added output purer than what was already there.

## Rule 4, and the correction I nearly shipped

My first sharing audit reported the sharing buckets as *purer* (alone 87%,
pairs 93%, groups 100%) and I was one commit from recording "the failure mode
was falsified in the opposite direction". Codex found three defects:

* it counted query-side `=` labels, so several added models matching the same
  reference each scored as a win;
* its predicate required one whole chain to contain the other, missing the
  ordinary case — `[A,C]` and `[B,C]` are both attestable by a read carrying
  `[C]` while neither contains the other;
* its grouping joined a model to the first group it matched and never merged
  the rest, so `[A,C]`, `[B,C]`, `[C]` came out as a pair plus a singleton.

Rebuilt with distinct-reference counting, any-shared-suffix linking and
union-find, "alone" drops from 178 models to 107 and the picture reverses:

| sharing bucket | models | new truth | marginal Pr |
| --- | ---: | ---: | ---: |
| alone (1) | 107 | 96 | **90%** |
| pair (2) | 64 | 48 | **75%** |
| group (≥3) | 61 | 50 | **82%** |

(replicated on `c_trunc5p50`: 91 / 79 / 82.)

So the shared-suffix buckets yield **8–15 pp less**, consistent with the
predicted vulnerability rather than contradicting it — and still clearing
break-even by 3.5–4×. The honest statement is that the vulnerable regime did
not cost us *here*, not that the vulnerability is absent. Two limits on even
that: the audit sees only models that survived the whole selection cascade, so
the GTF floor may already be removing the weak members of ambiguous families
(a reason to leave the floor untouched in the same change), and it groups by
structural ambiguity rather than observed shared reads.

## The configuration trap this uncovered

Adding the key to the profile did nothing. `resolve_profile_values` returns a
dict and the CLI must lift each key into the local it passes to
`PipelineConfig`; the profile entry existed, 3 was computed, and click's flag
default of 0 is what ran. Nothing raised. The only symptom was a verification
run coming out byte-identical to the arm it was supposed to differ from.

The guard is now behavioural: ten tests through the real click command
asserting on the config the runner is handed, including one that iterates every
profile and checks every key arrives. Verified to fail when the read-back line
is deleted. A textual source check — which is what I wrote first — passes on a
key that appears in a comment, in dead code, or in an assignment that is later
overwritten, and would not have caught this.

## Next lever

A `--selection-audit` pass on the human `full` cell under the new default is
running, to attribute the *remaining* loss now that the gate is no longer the
owner. Per Codex, the accounting to build from it is
`exempted → reaches global selection → GTF floor drop → final distinct truth`,
splitting floor casualties into zero-mass and (0,1). The floor stays fixed in
this change: it may be part of why the exemption is safe.
