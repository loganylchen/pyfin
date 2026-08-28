# PyFIN Repository Map

Snapshot: branch `dev`, clean baseline commit `97029eae` (`feat(pipeline): add
validated profiles and family-aware selection`), plus the v9 abundance-refit
working tree with source SHA-256 `42e1210d...7954`, updated through 2026-08-27.
The baseline branch is synchronized with `origin/dev`; this document describes
the current working tree, not only the last pushed commit.

## 1. Scope and evidence

- Core index: 181 files from `fin/`, `tests/unit/`, `tests/integration/`,
  `benchmarks/`, README, packaging metadata, and current `.omc` plans/specs.
- Experiment/provenance index: 262 tracked or pending Python, shell, SBATCH,
  and Markdown files under `experiments/`. Generated result trees and raw
  signal/alignment data were deliberately not content-scanned.
- Infrastructure index: 7 CI, Docker, history, and repository-policy files.
- Static architecture map: 67 Python modules, 115 internal import edges, 83
  formal unit/integration test files, and direct test-import relationships.
- The repository requests the `code-review-graph` MCP, but those tools were not
  exposed in this session. The fallback was a Git-bounded inventory plus AST
  import/test analysis. `.code-review-graph/` itself occupies about 28 MB.
- Source at `HEAD` is treated as authoritative where README, `record.md`,
  `PRODUCTION_STATE.md`, plans, comments, or old experiments disagree.

The pushed baseline has 551 tracked files; this v9 batch adds eight source,
test, audit, and SBATCH files. Generated experiment data was removed from Git
history before the prior push. Roughly 25 GB of reference genomes, reads,
comparison output, and the 1.87 GB production SIF remain available locally as
untracked artifacts; their generating scripts and roles are mapped without
committing or scanning every byte.

## 2. Product boundary

PyFIN (`py-fin` 0.1.0, alpha, Python >=3.8) is a nanopore Direct RNA-seq
transcriptome assembler. It combines sequence alignment, raw-signal evidence,
structural candidate generation, and probabilistic assignment to emit known,
novel, and optionally fusion transcripts.

Primary inputs:

- Spliced read alignments: BAM/SAM through `pysam`.
- Optional annotation: GTF/GFF.
- Genome sequence: FASTA.
- Basecalled reads: FASTQ.
- Raw signal: SLOW5/BLOW5 or POD5; FAST5 readers remain available as library I/O.

Primary outputs:

- `assembly.gtf`: surviving transcripts with abundance-related attributes.
- `scores.tsv`: candidate identity, locus, source, abundance, confidence,
  coherence/discrimination/combined scores, hard-read count, TPM, and `max_R`.
- `fusions.bedpe`: fusion breakpoints when fusion detection is enabled.
- `abundance_refit.json`: final-survivor mass conservation, orphaned reads,
  and abundance/TPM shift diagnostics for named profiles.
- Optional unfiltered diagnostic TSV and persisted responsibility matrices.

Internal genomic intervals, exons, and introns are represented as 0-based,
half-open coordinates. GTF is read/written as 1-based inclusive and converted at
the I/O boundary. Strand is always explicit where biological interpretation
depends on transcript orientation.

## 3. Runtime architecture

The current execution path is:

```text
fin CLI
  -> PipelineConfig construction
  -> PipelineRunner.setup()
       GTF reader + genome FASTA + signal reader
  -> PipelineRunner.run()
       isolated strand-aware genomic intervals
       -> process_interval() serially or through spawn workers
          -> candidate discovery (GTF + novel)
          -> optional fusion discovery and merge
          -> pre-assignment canonical and junction-dominance gates
          -> one of argmax / m1_em / m2_em / cluster assignment
          -> m2 interval-level structural selection
          -> per-interval QuantResult objects + sparse responsibility ledger
       -> commutative cross-interval aggregation
       -> optional unfiltered diagnostic
       -> global selection cascade
       -> finalized junction consensus/merge
       -> final-survivor responsibility refit + orphaned-mass audit
       -> gene-id resolution and GTF/TSV/BEDPE writers
  -> PipelineRunner.cleanup()
```

### 3.1 Interval construction

`fin.io.interval_manager` builds non-overlapping, strand-separated
`GenomicInterval` objects from BAM reads and optional GTF records. It tracks
read counts but not gene/transcript IDs. Reads identified as chimeric/fusion
reads are excluded from ordinary interval generation. Read extraction and GTF
annotation extraction are interval scoped.

### 3.2 Candidate construction

The production generation default is `chain_cluster_discovery=True` with
`clustering="families"`:

1. Extract read intron chains and transcript 3-prime positions.
2. Group reads structurally, with 6 bp junction wobble and a 70 bp cassette
   exon allowance.
3. Build `ChainFamily` groups by single-linkage structural relations.
4. Attach zero-read GTF hypotheses without allowing them to bridge families.
5. Apply explicit exact-subchain `collapse()` with a read-span guard.
6. Defer mono-exon resolution until after EM by default.
7. Optionally use canonical splice alternatives, de-novo junction graphs, or
   family-local path exploration. Some of these remain experimental/off by
   default.

`TranscriptCandidate.source` is one of `gtf`, `novel`, or `fusion`. Fusion
candidates carry two explicit breakpoints and are protected from ordinary
candidate collapse. `CandidateSet` also carries read sequences, family/cluster
membership, shadow-chain provenance, and observed read spans when available.

`fin.scoring.candidate_align` and `fin.candidates.diff_regions` are the two
newest foundations. They have tests and experiment consumers, but no inbound
import from the production `fin` graph at this snapshot. The same is true of
`fin.candidates.explore`. They should not be described as active production
stages until wiring is added.

### 3.3 Fusion path

Fusion detection is off by default. When enabled:

1. `fusion.chimeric` collects SA-tag/chimeric reads and re-aligns clipped arms
   with a cached genome-wide mappy aligner.
2. `fusion.arm_assembly` clusters breakpoint-compatible reads and infers each
   arm's splice structure.
3. `fusion.stitch` joins arm variants across the breakpoint into candidates.
4. `fusion.detect` orchestrates F1 -> F2 -> F3.
5. Candidates merge into the ordinary `CandidateSet`, retain breakpoint
   metadata through quantification, and are written to BEDPE.

The current soft-clip threshold for identifying a fusion read is 250 bp.

### 3.4 Evidence and scoring layers

- M1 sequence evidence: per-read/per-candidate minimap2/mappy alignment score,
  converted to a distance. Large-indel rejection and preset selection are
  centralized in `mappy_score.py` and `mappy_preset.py`.
- M1 redesign foundation: `candidate_align.py` calculates family-scoped
  compatibility/goodness, but is not yet production-wired.
- M2 signal evidence: krill eventalign over candidate-discriminating internal
  junction windows. In the default `m2_em` path, M2 only refines reads tied for
  best M1 alignment when the signal window is scorable.
- M3 read-to-read coherence was removed from production; its prototype and
  focused tests are preserved under `experiments/m3_coherence/`.
- M4 whole-read/diff-region DTW remains experimental with no production wiring.
  The structural `diff_regions.py` engine and older `diff_region_dtw.py`
  scoring implementation are distinct components.
- PolyA evidence: whole-read krill polyA estimation plus 5-prime proximity.

`fin._dtw` currently exposes a krill-backed compatibility API. It does not
import the `_cuda_dtw` extension that `setup.py` still builds. Host-side
`is_available()` was false; the intended production environment supplies krill
inside the SIF.

### 3.5 Assignment and quantification

Four `quant_mode` values exist:

- `argmax`: hard/split M1 assignment without signal or EM.
- `m1_em`: EM from M1 mappy distance, with the coherence term disabled.
- `m2_em` (default): one mappy AS pass determines kept reads and ties; junction
  NLL refines only tied cells; EM produces soft responsibilities; hard
  assignments and quantification are derived afterward.
- `cluster`: assignment is scoped inside candidate families/clusters, with an
  optional M2 leg.

The generic EM implementation supports `sigma`, `beta`, priors, optional
abundance feedback, optional effective-length normalization, CPU numpy, and GPU
CuPy. `Assigner` deliberately performs assignment only; candidate dropping is
owned by `selection.py`. Its immutable `QuantOutput` is the boundary between
assignment and interval selection.

`QuantResult` carries abundance, confidence, hard-read count/IDs, source,
locus/exons, gene ID, coherence/discrimination/combined scores, full-length
fraction, and optional fusion breakpoints. Cross-interval aggregation is
commutative and preserves these fields.

### 3.6 Selection ownership and order

Pre-assignment selection:

1. Canonical gate (default on): remove novel multi-exon candidates containing
   non-canonical junctions. Accepted motifs default to GT-AG, GC-AG, AT-AC.
2. Junction-dominance gate (default off): remove novel candidates whose
   junction evidence is weak or displaced by a nearby dominant junction.

`m2_em` interval selection can apply, in order, novel/guided junction support,
M2 support, containment-cluster, containment-shadow, deferred mono-exon,
structural-wobble, and related GTF-guarded rules. Dropped candidates are set
filtered; old precision-lever plans correctly note that most drops do not
redistribute their soft mass unless an explicit fold operation does so.

Global post-aggregation selection applies source-aware abundance floors,
family-aware isoform fraction (overlap fallback), soft-mass/hard-read ratio, optional mono-exon filtering,
full-length fraction, and optional polyA/5-prime support. Fusion candidates and
GTF candidates have explicit exemptions in several filters. Real profiles then
apply finalized-model read-supported junction consensus and mass-preserving
structural merges; SIRV/custom leave it off by default.

`finalize.py` does not select candidates. It resolves gene IDs and delegates to
the standalone GTF, TSV, and BEDPE writers. This separation is intentional so
output formatting cannot change the survivor set.

### 3.7 Parallelism

Interval work is process-parallel with multiprocessing `spawn`:

- `threads` is the total worker count.
- `gpu_workers` statically partitions workers into GPU and CPU runners, bounding
  the number of CUDA contexts.
- Each worker creates one runner in its initializer and reuses it.
- Only deterministic interval indices cross the process boundary.
- CUDA, pysam, POD5, and inherited BLAS pools are the reasons `fork` is avoided.
- Krill may fall back from GPU to CPU per worker after initialization/OOM.

## 4. Core module catalog

`Direct tests` counts formal unit/integration files that import the module
directly. It is not a coverage percentage: runner-level and mocked transitive
tests can exercise modules with a zero in this column.

| Module | LOC | Direct tests | Responsibility |
| --- | ---: | ---: | --- |
| `fin` | 55 | 2 | Package initialization, version, GPU capability logging |
| `fin.__main__` | 6 | 0 | `python -m fin` CLI bridge |
| `fin._dtw` | 180 | 1 | Krill-backed DTW compatibility API |
| `fin.ablation` | 1 | 0 | Ablation namespace |
| `fin.ablation.io` | 131 | 1 | Ablation summary/per-row TSV writers |
| `fin.ablation.mappy_argmax` | 230 | 0 | M1-only assignment and responsibilities |
| `fin.ablation.runner` | 126 | 1 | Active M1/M2 quantification ablation rows |
| `fin.analysis` | 14 | 0 | Analysis namespace |
| `fin.analysis.abundance_refit` | 376 | 1 | Final-survivor responsibility refit and mass diagnostics |
| `fin.analysis.assignments` | 285 | 4 | Coherence-aware EM |
| `fin.analysis.clustering` | 305 | 1 | Legacy 3-prime BAM clustering workflow |
| `fin.analysis.quantification` | 652 | 17 | Quantification, family-aware fraction, TPM, filters, aggregation |
| `fin.candidates` | 21 | 0 | Candidate namespace |
| `fin.candidates.canonical` | 131 | 2 | Motif parsing and canonical splice alternatives |
| `fin.candidates.chain_cluster` | 631 | 4 | Read-chain/family clustering and collapse |
| `fin.candidates.dataclasses` | 123 | 32 | IntronChain, family-tagged TranscriptCandidate, CandidateSet |
| `fin.candidates.diff_regions` | 322 | 1 | Cluster structural difference regions |
| `fin.candidates.discovery` | 935 | 8 | Per-interval discovery and stable candidate/family IDs |
| `fin.candidates.explore` | 167 | 1 | Family-local graph path enumeration |
| `fin.candidates.intron_chains` | 582 | 3 | Chain extraction, snapping, grouping, representatives |
| `fin.candidates.isoform_recovery` | 120 | 1 | Post-EM 5-prime/TSS peak recovery |
| `fin.candidates.junction_graph` | 193 | 1 | De-novo intron graph and chain assembly |
| `fin.cli` | 583 | 1 | Click CLI, named profiles, and manifest config mapping |
| `fin.fusion` | 46 | 0 | Fusion API exports |
| `fin.fusion.arm_assembly` | 285 | 2 | Fusion arm clustering and splice inference |
| `fin.fusion.chimeric` | 255 | 3 | Chimeric reads and clipped-arm re-alignment |
| `fin.fusion.detect` | 71 | 0 | Fusion stage orchestrator |
| `fin.fusion.stitch` | 148 | 1 | Cross-breakpoint candidate stitching |
| `fin.io` | 112 | 0 | Format API and nominal lazy loaders |
| `fin.io.interval_manager` | 655 | 13 | Strand-aware intervals and interval extraction |
| `fin.io.io_bam` | 407 | 0 | BAM/SAM parser |
| `fin.io.io_bed` | 384 | 0 | BED parser |
| `fin.io.io_bedpe` | 75 | 2 | Fusion BEDPE writer |
| `fin.io.io_fast5` | 221 | 0 | FAST5 reader |
| `fin.io.io_fasta` | 333 | 2 | FASTA parser |
| `fin.io.io_fastq` | 201 | 0 | FASTQ parser |
| `fin.io.io_gtf` | 651 | 2 | GTF/GFF reader and writer |
| `fin.io.io_pod5` | 363 | 0 | POD5 reader |
| `fin.io.io_read_manager` | 675 | 0 | Integrated BAM/GTF/FASTA read subsets |
| `fin.io.io_slow5` | 303 | 0 | SLOW5/BLOW5 reader and pA conversion |
| `fin.io.io_tsv` | 58 | 3 | Candidate scoring TSV writer |
| `fin.pipeline` | 9 | 2 | Pipeline namespace |
| `fin.pipeline.assignment` | 592 | 1 | M2 assignment, tie NLL, and per-read batch abstention |
| `fin.pipeline.cluster_quant` | 212 | 1 | Within-family assignment |
| `fin.pipeline.config` | 765 | 14 | 127-field configuration and named-profile contract |
| `fin.pipeline.evidence` | 93 | 0 | Reusable per-interval junction evidence |
| `fin.pipeline.finalize` | 109 | 0 | Diagnostics, gene IDs, output wiring |
| `fin.pipeline.junction_snap` | 235 | 1 | Finalized junction consensus, merges, and ledger redirects |
| `fin.pipeline.parallel` | 139 | 0 | Spawn worker lifecycle, partitioning, and ledger transport |
| `fin.pipeline.runner` | 1396 | 10 | End-to-end interval/pipeline/refit orchestrator |
| `fin.pipeline.selection` | 619 | 0 | All pre/interval/global survivor decisions |
| `fin.scoring` | 16 | 0 | Scoring namespace |
| `fin.scoring.candidate_align` | 182 | 1 | Family-scoped M1 goodness foundation |
| `fin.scoring.diff_region_dtw` | 821 | 1 | Difference-region M4/coherence scoring |
| `fin.scoring.em_inputs` | 126 | 1 | Coherence-free M1/M2 matrix subset construction |
| `fin.scoring.eventalign_parser` | 293 | 3 | Legacy f5c eventalign TSV parsing |
| `fin.scoring.krill_aligner` | 119 | 1 | Shared krill backend and CPU fallback |
| `fin.scoring.krill_tiebreak` | 321 | 0 | In-memory signal tie-break |
| `fin.scoring.m2_junction_nll` | 1464 | 12 | Junction windows, NLL, structural/support drops |
| `fin.scoring.mappy_distance` | 102 | 1 | M1 read-candidate distance matrix |
| `fin.scoring.mappy_preset` | 20 | 0 | Central mappy preset |
| `fin.scoring.mappy_score` | 56 | 0 | Alignment score reconstruction/indel rejection |
| `fin.scoring.polya` | 99 | 1 | Whole-read krill polyA estimation |
| `fin.scoring.signal_dtw` | 300 | 1 | Signal slicing and pairwise DTW |
| `fin.utils` | 7 | 0 | Utility namespace |
| `fin.utils.log_config` | 131 | 1 | File/console logger configuration |
| `fin.utils.sequences` | 28 | 1 | Reverse complement |

The largest and most change-sensitive modules are
`scoring/m2_junction_nll.py` (1,464 LOC), `pipeline/runner.py` (1,396),
`candidates/discovery.py` (935), `scoring/diff_region_dtw.py` (821), and
`io/io_read_manager.py` (675). Responsibility extraction reduced runner size,
but these remain the primary review hotspots.

## 5. Configuration contract

`PipelineConfig` has 127 fields. The Click CLI exposes 102 options. Programmatic
defaults intentionally differ from raw Click and named-profile defaults: the
dataclass uses `min_abundance=0.0` and `min_gtf_abundance=0.0`, raw Click uses
3.0 and 1.0, and the default `real-drna` profile resolves strict >1.0 plus its
validated mono/junction policies. Programmatic direct construction does not
inherit named-profile settings unless it opts in.

Complete field/default groups at this snapshot:

```text
Inputs/work:
  bam_path=<required>; gtf_path=None; genome_fasta_path=''; fastq_path='';
  signal_path=''; signal_format='slow5'; work_dir='./pyfin_work'

Generation:
  three_prime_threshold=24; max_gap=0; min_novel_reads=1;
  chain_cluster_discovery=True; clustering='families';
  chain_cluster_wobble_bp=6; chain_cluster_cassette_max_exon_bp=70;
  chain_cluster_fold_monoexon=True; chain_cluster_fold_span_guard=True;
  mono_resolve_post_em=True; mono_resolve_min_reads=2;
  mono_resolve_slop_bp=10; canonical_search_bp=4;
  max_chains_per_read=16; canonical_gate=True;
  canonical_motifs=('GT-AG','GC-AG','AT-AC')

Global selection:
  min_abundance=0.0; floor_gtf_abundance=False; min_gtf_abundance=0.0;
  min_isoform_fraction=0.01; isoform_fraction_locus='family';
  post_selection_refit=False (named profiles=True; effective resolved by validate);
  max_soft_mass_ratio=2.0;
  min_fulllen_fraction=0.1; fulllen_window_bp=25; fulllen_min_reads=4;
  min_polya5p_reads=1; polya5p_window_bp=25; min_polya_length=10.0;
  polya5p_exempt_gtf=True; drop_mono_exon_novel=False;
  min_mono_exon_reads=0; min_mono_exon_length=0;
  junction_snap=False; junction_snap_tolerance=6;
  junction_snap_min_support=2; junction_snap_min_ratio=2.0

EM/quantification:
  em_sigma=1.0; em_max_iter=1000; em_tol=1e-4;
  abundance_feedback=False; abundance_length_norm=False;
  quant_mode='m2_em'; cluster_llr_threshold=2.0;
  cluster_min_support=1.0; cluster_m1_tie_margin=20.0;
  cluster_use_m2=True; em_max_iter_override=None;
  enable_score_filter=True; tiebreak_ambig_threshold=0.9;
  krill_tiebreak=False; krill_pore='rna002'

M2 and structural selection:
  m2_tiebreak=True; m2_tiebreak_junction_k=10; m2_tiebreak_margin=1e-9;
  m2_diff_cover_gate=True; m2_diff_cover_margin=0.5;
  m2_cluster_recheck=True; m2_cluster_recheck_bp=20;
  m2_cluster_recheck_fraction=0.15;
  m2_cluster_recheck_cassette_max_exon_bp=70;
  m2_cluster_recheck_novel_displaces_gtf=True;
  m2_cluster_recheck_gtf_min_jct_reads=1;
  m2_cluster_recheck_jct_tol=0; containment_collapse=False;
  containment_3p_tol_bp=20; containment_min_abundance_ratio=1.0;
  containment_cluster=True; containment_cluster_wobble_bp=6;
  containment_cluster_min_ab_ratio=0.3;
  containment_cluster_min_read_ratio=0.3;
  containment_cluster_max_shadow_reads=10

Junction/de-novo gates:
  novel_junction_min_reads=2; novel_junction_reads_tol=2;
  guided_junction_min_reads=0; guided_junction_reads_tol=2;
  denovo_wobble_tol=0; denovo_wobble_shadow_ratio=0.5;
  denovo_graph=False; denovo_graph_tol=6; denovo_graph_min_edge_reads=2;
  denovo_graph_tss_brake=True; denovo_graph_tss_tol=20;
  denovo_graph_tss_min_reads=3; denovo_graph_tss_frac=0.4;
  junction_dominance_filter=False; junction_dominance_min_reads=2;
  junction_dominance_window_bp=20; junction_dominance_tol_bp=2;
  m2_support_gate=True; m2_support_gate_tie=True;
  m2_tie_scoregate_split=True

Scoring/runtime/output:
  score_alpha=0.5; prior_weight_cap=10.0; use_prior=True; use_gpu=True;
  max_reads_per_interval_for_dtw=2000; signal_normalize=True;
  threads=1; gpu_workers=0; output_gtf=None; output_tsv=None;
  output_bedpe=None; max_reads=None; persist_R_matrix=True;
  write_unfiltered_scores=False

Fusion:
  fusion_enabled=False; fusion_min_support=2; fusion_max_dist=500;
  fusion_flank_bp=500; fusion_max_internal_gap_bp=30
```

CLI overrides worth treating as compatibility-sensitive:

- `min_abundance=3.0`, `min_gtf_abundance=1.0`.
- `use_gpu=True`, `quant_mode='m2_em'`, `clustering='families'`.
- Canonical, M2 tie-break, diff-cover, cluster-recheck, M2 support,
  containment-cluster, mono-resolve, span-guard, and responsibility persistence
  default on.
- M3 coherence, de-novo graph, containment-collapse, junction dominance,
  mono-exon drop, abundance feedback, and length normalization default off.
- Full-length and polyA/5-prime filters have nonzero CLI defaults, but comments
  in production notes say they are SIRV-oriented and may be disabled for real
  direct-RNA operation. Treat source/config and run manifests as the truth for a
  given experiment.

## 6. Packaging, containers, and CI

Declared core libraries include numpy, pandas, scipy, pysam, mappy,
ont-fast5-api, pyslow5, pod5, h5py, Click, tqdm, matplotlib, seaborn, and YAML.
CuPy/numba are optional GPU extras. Krill is operationally required for current
signal/DTW paths but is supplied by container images rather than declared in
`pyproject.toml`.

Known metadata drift:

- `pyproject.toml` requires `pyslow5>=0.3.0`; `setup.py` requires >=1.0.0 and
  explains that older versions silently lose reads during pA conversion.
- `requirements.txt` omits `mappy`, while pyproject/setup include it.
- Pyproject installs the `fin` command; setup.py installs both `fin` and
  `pyfin`.
- `setup.py` says it builds f5c integration and builds `_cuda_dtw`; the current
  Python DTW facade uses krill and never imports `_cuda_dtw`.
- `fin.io` advertises lazy optional loaders, but importing the package can still
  reach `interval_manager -> io_bam` and require pysam eagerly.

Runtime used by production/Slurm scripts:

```bash
SIF=experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
singularity exec --nv -B /SSD "$SIF" \
  env PYTHONPATH=/SSD/logan/dev/pyfin \
  /usr/bin/python3.10 -m fin.cli ...
```

There is no repository-local venv. Host Miniforge has numpy/pysam but lacks
pytest/mappy; the host `pytest` executable uses `/usr/bin/python3`, which lacks
numpy/pysam. Host test results are therefore meaningless unless the interpreter
is controlled explicitly.

CI/infrastructure:

- `python-package.yml` tests Linux/macOS with Python 3.8-3.11 and runs lint,
  packaging, and coverage. It triggers on `main` and `develop`, while the active
  branch is named `dev`; ordinary pushes to `dev` do not trigger that workflow.
- `docker-build.yaml` runs on every branch, builds/pushes GPU and CPU images,
  and performs import/backend smoke tests.
- CPU, GPU, and krill-specific Dockerfiles exist. The production SIF is tracked
  under experiments rather than rebuilt locally during normal development.

## 7. Test map and observed health

Pytest configuration formally targets `tests/unit` and `tests/integration` with
strict markers. Numerous scripts directly under `tests/` are diagnostic,
benchmark, visualization, or manual real-data programs and are not part of the
configured formal suite.

Observed in the repository SIF with the final profile-integration source:

- Unit: 1032 passed, 1 skipped in 9.76 s (`singularity --nv`).
- Integration: 15 passed, 3 skipped in 1.95 s.
- Retired M3 prototype: 5 passed in 1.05 s.
- CLI parameter tests now create their own existence-only temporary inputs, so
  they reach `threads`/`gpu_workers` validation instead of failing on absent
  ignored fixtures.
- `PipelineConfig.validate()` is called after profile resolution and before
  runner setup, with failures translated to Click usage errors.
- Benchmark smoke tests and `run_benchmark.sh` use `sys.executable` or a
  discovered Python 3 executable rather than assuming bare `python`.
- Skips are the explicitly optional large-fixture/GPU-equivalence cases whose
  assets are not part of the formal test checkout. Fresh SIRV and full real p00
  end-to-end runs provide direct pipeline evidence beyond those skips.

## 8. Branch evolution and active work

The 19 commits ahead of `origin/dev` tell one coherent architecture story:

1. `c49fd9a1`: generation-side intron-chain clustering.
2. `828565a7`: fusion soft-clip threshold raised to 250 bp.
3. `7420fe6a`: production chain-cluster body and evidence-layer foundation.
4. `bc39ce8b`: extract M2 tie NLL to assignment.
5. `b01bcf61`: extract M2 assignment into `Assigner`/`QuantOutput`.
6. `6d99c5b6`: move interval selection into `selection.py`.
7. `e1698ae9`: move global selection into `selection.py`.
8. `cda1a5b4`: move pre-assignment gates into `selection.py`.
9. `1f754007`: move output finalization/writers into `finalize.py`.
10. `84241b10`: add grouping-only `cluster_families`.
11. `d277b4ee`: attach non-bridging, zero-read GTF hypotheses.
12. `65510187`: add family-local intron-graph exploration.
13. `614f985b`: add explicit exact-subchain family collapse.
14. `61c2f650`: wire family/collapse behind the clustering option.
15. `e8c16417`: SIRV clustering validation harness.
16. `35d3ec99`: span guard makes family collapse reproduce read-chain generation.
17. `48343336`: switch production default to families plus mono resolve.
18. `871f08d4`: add family-scoped read-candidate goodness foundation.
19. `a2c4be4e`: add cluster structural difference regions.

This range also tracks very large validation artifacts, genomes, reads, a SIF,
gffcompare output, and review transcripts. The source-level architecture change
is much smaller than the multi-gigabyte Git range suggests.

Current working tree before this map was created:

- One tracked user modification: `.gitignore` (+2 lines).
- 230 untracked files.
- Tool/config files for code-review-graph across multiple coding clients.
- GENCODE/SIRV canonical and gate-validation output/manifests.
- A large `experiments/wobble_heya8/` investigation into family differences,
  anchored two-hypothesis signal NLL, mappy goodness, read support, AUROC,
  cluster-bootstrap confidence intervals, and visual diagnostics.

The current `.omc` precision-lever plan is marked pending approval and asks for
default-off, novel-only containment mass-fold, perfect-junction, and mono-read
support levers. Parts overlap capabilities already present in history. Its most
important verified constraint remains valid: `io_bam.alignment_to_dict` does
not expose NM/MD/tag access, so a truly mismatch-free "perfect alignment" gate
needs new data plumbing.

`experiments/prod_validation/PRODUCTION_STATE.md` is dated 2026-07-22 and
predates most of the architecture extraction and the latest scoring/diff-region
commits. `record.md` contains valuable experiment decisions, but statements
such as M2 primary assignment being a dead end conflict with the current
`quant_mode='m2_em'` default. Both are historical evidence, not current source
truth.

## 9. Known gaps and review hotspots

1. Packaging and runtime declarations have drifted around krill, pyslow5,
   mappy, entry points, and the unused CUDA extension build.
2. Historical product truth remains distributed across source and experiments;
   `PROFILE_OPTIMIZATION.md` plus source-hashed manifests now define the current
   profile evidence, but older mutable notes still need archival labels.
3. `experiments/` still stores 25+ GB locally, but generated data is no longer
   tracked. Repository hygiene now depends on adding durable ignore rules so a
   future broad `git add` cannot reintroduce those artifacts.
4. `candidate_align`, `diff_regions`, and `explore` are tested foundations but
   not production-wired; they have no winning end-to-end evidence.
5. Several I/O modules and extracted pipeline layers have no direct-import test
   file. Transitive tests exist, but direct contracts are less explicit.
6. The two largest logic modules, `m2_junction_nll.py` and `runner.py`, still
   concentrate many policy interactions and deserve impact analysis for every
   change.
7. CI branch names still do not include `dev`; benchmark Python portability is
   fixed, but ordinary pushes to the active branch can miss package CI.
8. Same-chain endpoint isoforms and family-level novel gene IDs remain open.
   Production beta=0 post-drop abundance refitting and stable IDs are landed.
9. Real profile optimization now includes the tuning sample plus two independent
   H9 samples. Balanced/precision overrides and manifests remain necessary
   because no finite benchmark proves universal transcriptome optimality.

## 10. Change guidance

For future work, use this ownership sequence:

- Candidate-set changes: `candidates/*`, then discovery/family tests, then
  runner integration and real-data candidate-count/recall ablations.
- Assignment changes: `pipeline/assignment.py` and `scoring/em_inputs.py`; do
  not drop candidates there.
- Survivor-policy changes: `pipeline/selection.py` and the corresponding
  structural helper in `m2_junction_nll.py`/quantification.
- Finalized structural correction: `pipeline/junction_snap.py`; preserve read
  IDs and emit absorbed-to-representative redirects whenever models merge.
- Final abundance changes: `analysis/abundance_refit.py`; keep existence frozen,
  conserve one unit per assignable read, and report selection-orphaned mass.
- Output changes: `pipeline/finalize.py` plus standalone writer tests; do not
  change selection there.
- Parallel changes: preserve spawn, deterministic interval ordering, bounded
  CUDA contexts, and commutative aggregation.
- Signal changes: verify both missing-signal/CPU fallback and krill GPU paths;
  report the actual backend, not the requested backend.
- Any default change: compare dataclass, CLI, ablation builders, Docker/SIF run
  manifests, README, and production notes separately.

This map is a navigation artifact. Exact function bodies and experiment scripts
remain indexed under `pyfin-core-repository`, `pyfin-experiment-code`, and
`pyfin-infra-and-history` for symbol-level retrieval.
