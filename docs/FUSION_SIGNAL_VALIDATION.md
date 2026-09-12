# Fusion junction signal and final-assembly evidence

This extends the earlier sequence-only implementation. It changes only the
independent fusion modules and post-assembly CLI orchestration. Ordinary
assembly discovery, assignment, selection, EM and refit remain unchanged.

## Evidence contract

1. Load the **completed** `assembly.gtf` and `junction_support.tsv`, read-only.
   A junction must occur in a surviving transcript and have at least three
   indel-clean reads. Overlapping interval rows use the maximum count, not a
   sum. Negative-strand coordinates are converted to transcript orientation.
2. These junctions nominate nearby alternatives and annotate structural
   support. Aggregated counts cannot exclude the current fusion read; they
   are **not independent signal confirmation**. A supported donor and acceptor
   on separate genes do not constitute fusion-edge support. No assembly count
   is added to complete-read fusion support.
3. Sequence correction and whole-transcript validation still run first.
   Assembly nominations have the explicit `final_assembly_clean` origin.
4. Signal validation uses krill eventalign and the existing assembly helpers
   `diff_junction_windows`, `_zrecords`, `tx2genome_array`, and
   `_mean_nll_in_window(..., reduce="sum")`. Alignment slicing uses the shared
   `score_hit` reconstruction and per-variant transcript start offsets.
   It does not use raw-DTW SSE or the experimental neighbor-core scorer.
5. Each ordered chromosome/strand arm receives a disjoint forward virtual
   interval. This preserves exon lengths and sequence order for internal
   introns and N-arm fusion edges. Virtual coordinates are never emitted as
   genomic GTF coordinates. Signal audit windows are explicitly marked virtual.
6. At each junction, compare the selected complete template against anchored
   coordinate alternatives, holding the other junctions fixed. This reuses
   assembly's scoring method, **not its exact-AS-tie nomination gate**.
   Every hypothesis must return a successful record, finite scores, spanning
   coverage and at least three finite events on each shared flank.
7. A decisive conflict removes that already-enumerated path alternative.
   Another sequence-validated alternative may survive, resolving coordinates
   without discarding the molecule. If none survives, the molecule supplies no
   complete-chain vote. Signal never promotes a sequence-invalid path or
   assembles a new combination of independently winning edges.
8. Missing/failed signal, duplicate/missing hypothesis records, insufficient
   coverage, nomination/hypothesis limits and weak margins explicitly abstain.
   Identical full-template alternatives conservatively abstain rather than
   resolving microhomology using alignment-dependent event placement.
   Abstention retains labeled sequence evidence, not a signal pass.

## Configuration and outputs

`--fusion` automatically receives `cfg.signal_path`, `cfg.krill_pore`,
`cfg.use_gpu`, `cfg.m2_summed_llr_flank` and `cfg.m2_summed_llr_margin`
after ordinary assembly completes. The real-drna replay uses rna002, CPU,
8 bp windows and a 1.0 summed-NLL margin. Fusion's additional defaults are
three finite events per shared flank and 128 hypotheses per molecule.
The clean-assembly threshold is three reads. These are defaults, not a claim
of independently calibrated fusion thresholds.

Standalone usage:

```python
from fin.fusion.pipeline import run_fusion_analysis
from fin.fusion.signal import FusionSignalSettings

result = run_fusion_analysis(
    "reads.bam", "reference.fa", "output",
    gtf_path="annotation.gtf", sequence_path="reads.fastq",
    signal_path="reads.blow5",
    signal_settings=FusionSignalSettings(pore="rna002", flank=8, min_margin=1.0),
    assembly_gtf_path="completed_assembly/assembly.gtf",
    assembly_support_path="completed_assembly/junction_support.tsv",
)
```

Both assembly-evidence paths must be provided together. Omitting `signal_path`
keeps the standalone API sequence-only. Missing explicitly supplied files fail
publication; an unavailable backend or absent read signal is instead audited
as abstention. All supplied evidence files are included in input snapshots.

- `junction_signal.jsonl`: every attempted junction, coordinate alternatives,
  local sums, event counts, margin, support/conflict/abstention, and whether the
  parent alternative survived. Successful backend calls report effective GPU.
- `fusion_reads.jsonl`: sequence evaluations, selected surviving paths,
  signal decisions, and assembly structural annotations. Original alternative
  indices are retained; `path_coordinates` disambiguates them after filtering.
- `junction_corrections.jsonl`: the preceding **sequence-stage** decisions;
  their original `not_scored` field is not the later signal verdict.
- `fusions.tsv`: `signal_status` distinguishes `signal_supported`,
  `partially_signal_supported`, `sequence_only_signal_abstained` and
  `not_scored`. A coordinate result may be `signal_resolved` even if invariant
  internal junctions had no signal contrast. Existing minimum-support and
  readthrough classifications still apply.
- `fusion_manifest.json`: requested settings, actual aggregate scoring status,
  assembly-evidence count, file snapshots, output hashes and source hash.

## Verification

The 27 added tests cover virtual-coordinate projection on both strands and
multi-arm/internal junctions; exact shared summed-NLL arithmetic; missing,
duplicate and failed records; insufficient/invalid events; signal conflicts;
surviving alternatives; identical-template abstention; resource limits;
real mappy slicing and per-read backend reuse; clean final-assembly membership;
duplicate-count handling; and file-to-output provenance/abstention.
All 27 pass. The broader unit/integration run reports 2360 passed, 27 skipped
and two existing `TestTieNllKrillUnavailable` failures: the native Python
environment lacks krill while those assembly tests import it unconditionally.
Both failures were already reproduced against the archived HEAD baseline in
the preceding validation. Actual signal validation uses the CPU krill container.

```bash
python -m pytest tests/unit/test_fusion_signal.py \
  tests/unit/test_fusion_assembly_evidence.py \
  tests/integration/test_fusion_independent_pipeline.py::test_pipeline_reads_final_assembly_evidence_and_audits_signal_abstention -q
```

Actual CPU krill 0.1.6 eventalign was also exercised on the existing 575-record
chrIS Sequin subset of `SGNex_K562_directRNA_replicate4_run1`, with ordinary
real-drna assembly settings. The completed assembly contributed 53 eligible
junctions. On molecule `86dba9e7-765e-44ea-a173-6abf9852d633`, fusion-edge
alternatives `(6949422, 6942091)` and `(6949423, 6942091)` produced summed NLL
41.033192940669544 (20 events) and 62.16819541814077 (35 events): margin
21.135002477471225. The first sequence-validated path survived; the second was
removed. The internal introns had no alternative and explicitly abstained.

This remains a **single-read, insufficient-support, readthrough/splicing
candidate**, not a validated biological fusion or a fusion accuracy benchmark.
The example demonstrates the signal path executes and affects adjudication;
it does not establish that the winning coordinate is ground truth.

Baseline HEAD, fusion-off and fusion-on runs have identical assembly hashes:

| Artifact | SHA-256 |
| --- | --- |
| `assembly.gtf` | `2a393d30f83189d39d772c5d672c042bd83eb47d1e672d910a1143c2082c6a9f` |
| `scores.tsv` | `99d8fe803007ae488f30907a11ca872a9b44a563dc31bca55162350c54fa637b` |
| `abundance_refit.json` | `5e6b5465231195590f0f719eec614404c17f021b81fc99682f0495b378d444a0` |

Local replay inputs, commands and logs are under
`/tmp/pyfin-fusion-cli-blyaq7b0/sequin/`; baseline is `final_baseline`, the cold
off replay is `signal_final_off`, and signal-enabled replays include
`signal_verified_current`. These local temporary artifacts are not portable test
fixtures. No comparative sensitivity, precision, F1, runtime or memory
improvement is claimed without a separate labeled benchmark.
