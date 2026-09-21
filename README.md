# PyFIN - Signal-Based Reference-Based Transcriptome Assembly

A Python package for nanopore Direct RNA-seq **transcriptome assembly with fusion
detection and per-transcript signal-based scoring**. Reference-based discovery of
known and novel isoforms; independent multi-gene fusion analysis; and bounded krill eventalign
that refines exact M1 ties before read-to-transcript assignment, yielding
signal-aware abundance / TPM per candidate.

## Where to read what

| You want to | Go to |
|---|---|
| Install it (**including the signal backend -- it is easy to miss**) | [Installation](#installation) below |
| Run it | [Quick Start](#quick-start) below |
| Prepare inputs, diagnose a run, interpret the outputs | [`docs/USAGE.md`](docs/USAGE.md) |
| Know what the algorithm actually does | [`docs/CURRENT_ALGORITHM.md`](docs/CURRENT_ALGORITHM.md) |
| Find any other document | [`docs/README.md`](docs/README.md) -- the index, including which documents are historical and what is known to be missing |

Several documents in this repository read as authoritative because each was,
on the day it was written. `docs/README.md` says which one is current for
which question; `ALGORITHM.md` and `PIPELINE_AS_IMPLEMENTED.md` in particular
are dated snapshots and carry banners saying so.

## Pipeline Overview

```
BAM + GTF + POD5/SLOW5/BLOW5 + FASTQ + genome FASTA
                |
        Phase 1: Candidate Discovery
                  - GTF + novel isoforms
                v
        Phase 2: M1 alignment + M2 junction eventalign
                |
                v
        Phase 3: Read-to-transcript assignment (CuPy / numpy)
                  - signal-refined soft responsibilities
                |
                v
        Phase 4: Existence selection + junction consensus
                |
                v
        Phase 5: Survivor-set abundance refit + outputs
                  - assembly.gtf       (with abundance attributes)
                  - scores.tsv         (per-candidate metrics + TPM)
                  - abundance_refit.json (mass/TPM diagnostics)
                |
                v
        Optional independent fusion analysis (--fusion)
                  - complete multi-arm read paths + junction correction
                  - normal/paralog mapping competition
                  - fusion/ tables, FASTA, arm paths and audits
```

### Phase Details

1. **Candidate Discovery** - Reference-based: GTF transcripts + novel isoforms via intron-chain extraction + 3' end clustering
2. **Signal Scoring** - M1 supplies a mappy structural mask; krill M2 junction NLL refines exact ties inside shared genomic disagreement windows
3. **Assignment** - Probabilistic read-to-transcript responsibilities with production read-by-read coherence disabled (`beta=0`)
4. **Existence selection** - Structural/support gates fix the final transcript set; real profiles then apply read-supported junction consensus
5. **Quantification** - Named profiles renormalize each assignable read over the final survivors, report orphaned mass, and write coherent abundance to the GTF and both abundance and TPM to the TSV
6. **Independent Fusion Analysis** *(optional, `--fusion`)* - Scans BAM/SA evidence and realigns full molecules and clipped arms; validates complete gene chains and junction corrections without inserting candidates into assembly or changing its read assignments.

## Installation

### 1. PyFIN itself

```bash
git clone https://github.com/loganylchen/pyfin.git
cd pyfin
pip install -e .
```

That gives you the assembly pipeline. `mappy` (the minimap2 Python binding)
and `pysam` come in as dependencies -- PyFIN calls them in-process and never
shells out to the `minimap2` or `samtools` CLIs.

### 2. krill -- required for signal scoring

**The shipped profile requires krill and fails without it.** `real-drna`
resolves `m2_metric = "summed_llr"`, and on that path the assignment layer's
`tie_nll` does an unguarded `import krill` (`fin/pipeline/assignment.py:369`)
*before* it checks whether any read actually needs tie refinement. A run
without krill therefore raises `ModuleNotFoundError` as soon as an interval
reaches assignment with reads still retained after candidate remapping --
serial or parallel alike (`fin/pipeline/parallel.py:128` propagates worker
failures). It does not silently degrade to a sequence-only result. Install
krill before you start.

Two exceptions, for completeness: an explicit `--m2-metric off` bypasses the
import, and an interval with no candidates, no read sequences, or no candidate
alignment clearing the AS filter returns before reaching it
(`fin/pipeline/assignment.py:964`, `:1009`). Neither is the default.

```bash
pip install krill --no-deps \
    --index-url https://loganylchen.github.io/krill-dist/simple/
```

krill is a faithful port of the f5c engine and supplies BOTH the in-memory
eventalign and the DTW. It replaces the external `f5c` CLI and the CUDA DTW
extension PyFIN used to ship; neither is needed or used any more.

### 3. Optional: GPU acceleration for the EM

```bash
pip install -e ".[gpu]"      # CuPy, for the assignment/EM step
```

CuPy accelerates the read-to-transcript EM; krill does the signal. They are
separate installs and either can be present without the other -- but they are
**not** independent of `--no-gpu`: that flag is passed to both. A CUDA-built
krill uses the GPU when one is visible and falls back to CPU on init failure
(`fin/scoring/krill_aligner.py:66`), and `--no-gpu` asks it for CPU directly.

### Verify the installation

```bash
# 1. the package imports
python -c "import fin; print(f'PyFIN v{fin.__version__}')"

# 2. the SIGNAL BACKEND is importable. The default profile needs it and will
#    raise ModuleNotFoundError during assignment if it is absent
python -c "from fin._dtw import is_available; print(f'krill backend: {is_available()}')"

# 3. optional GPU EM
python -c "from fin.analysis.assignments import CUPY_AVAILABLE; print(f'CuPy: {CUPY_AVAILABLE}')"
```

Check 2 reports whether **krill is importable**, which is what signal scoring
needs. It is not a CUDA probe: a CPU-only krill build reports `True`, and that
is correct -- it is installed and usable. Whether krill will actually use a
GPU depends on the build (`krill.built_with_cuda`), on a visible device, and
on `--gpu/--no-gpu`. Note that an import check is not a signal-health check;
there is no end-to-end smoke test that proves signal scoring actually ran
(see `docs/README.md`, known gaps).

## Quick Start

### Operating profiles

The CLI resolves a named scientific operating point before constructing
`PipelineConfig`. The product default is balanced `--profile real-drna`; use
`--profile real-drna-precision` for a higher-F1 support floor,
`--profile sirv` for synthetic SIRV benchmarking, and `--profile custom` for raw
option defaults. Any explicitly supplied CLI option overrides its profile value.

| Profile | Intended data | Key resolved behavior |
| --- | --- | --- |
| `real-drna` | Biological nanopore dRNA, balanced | Tight summed LLR; strict >1 abundance; novel mono >=5 hard reads; read-supported junction consensus; final survivor-set abundance refit |
| `real-drna-precision` | Biological nanopore dRNA, precision-prioritized | Balanced real settings plus generation support >=2 and final survivor-set abundance refit |
| `sirv` | SIRV synthetic benchmark | M2 auto-routes to summed LLR (1,4) without a usable GTF and mean with a guide; inclusive >=3 abundance; final survivor-set abundance refit |
| `custom` | Manual/legacy experiments | No profile overlay or refit; pass `--isoform-fraction-locus overlap` to reproduce pre-family selection |

Every run writes `run_manifest.json` beside the GTF. It records the selected
profile, explicit overrides, resolved configuration, source SHA-256, Git commit
when available, and result-changing environment variables. Refit-enabled runs
also write `abundance_refit.json` with assignable/alignment-unassigned reads,
selection-orphaned mass, conservation error, and abundance/TPM shifts.

```bash
# Biological sample (default profile shown explicitly for reproducibility)
pyfin --profile real-drna ...

# Precision-prioritized biological sample
pyfin --profile real-drna-precision ...

# SIRV: auto selects SUM when unguided and mean when a usable GTF is supplied
pyfin --profile sirv ...

# Explicit A/B override remains available
pyfin --profile sirv --m2-metric mean ...
```

### CLI — assembly only

```bash
pyfin \
    --bam reads.bam \
    --gtf annotations.gtf \
    --genome genome.fa \
    --fastq reads.fq \
    --signal signal.blow5 \
    --signal-format slow5 \
    --output-dir pyfin_out/
```

Outputs:
- `pyfin_out/assembly.gtf` — assembled transcripts with abundance, confidence, hard-read count, source, and diagnostic score attributes
- `pyfin_out/scores.tsv` — per-candidate metrics and TPM table
- `pyfin_out/candidate_evidence.tsv` — optional (`--candidate-evidence`):
  per-survivor inference-time observable features (junction support, family
  share, end support, canonicality, containment geometry); read-only
  diagnostic feeding the experimental evidence ranker (`--ranking-mode
  filter`, off by default)
- `pyfin_out/abundance_refit.json` — survivor-refit mass conservation and abundance/TPM shift diagnostics (named profiles)
- `pyfin_out/endpoint_refine.json` — optional (`--endpoint-refine`): endpoint
  modes, pair support, poly(A) support, and read routes per split model
- `pyfin_out/catalog_completion.json` — optional (`--catalog-mode complete`
  or `adaptive`): catalog additions, adaptive decision, exact mass-transfer
  map, and strict before/after conservation audit. In adaptive mode, source
  filtering alone does not reconstruct pre-transfer abundance.
- `pyfin_out/m2_contrasts/*.jsonl` — optional
  (`--m2-contrast-stats-jsonl`): per-attempt M2 contrast records (event
  counts, mean/sum NLL, margin, coverage, abstention reason)
- `pyfin_out/run_manifest.json` — resolved profile, overrides, configuration, source hash/commit, and result-changing environment

Experimental opt-in flags (off in every profile; each requires its own
validation before promotion):

| Flag | Effect |
| --- | --- |
| `--candidate-evidence` | write the observable-feature table (read-only) |
| `--ranking-mode filter` | drop novel candidates below the frozen evidence-ranker operating point before snapping/refit |
| `--endpoint-refine` | split novel multi-exon models into supported TSS/TES states, with mandatory post-split requantification |
| `--catalog-mode complete\|adaptive` | catalog-structure benchmark adapter only: full-GTF gffcompare rewards zero-abundance entries, but expression-conditioned rescoring shows large precision loss; never enabled by a profile and not an expression-performance mode (the supporting study is an internal note, not part of this repository) |
| `--m2-metric sqrt_count_mean_llr` | event-count-rescaled tight-window M2 contrast |
| `--m2-contrast-stats-jsonl` | machine-readable M2 contrast records |
| `--no-lazy-genome` | restore the eager whole-genome dict (higher worker RSS) |

### CLI — assembly + fusion detection

```bash
fin \
    --bam reads.bam --gtf annotations.gtf --genome genome.fa \
    --fastq reads.fq --signal signal.blow5 \
    --output-dir pyfin_out/ \
    --fusion --min-support 2 --max-dist 6 --flank-bp 35
```

Adds results under `pyfin_out/fusion/`. Ordinary `assembly.gtf`, `scores.tsv`
and `abundance_refit.json` remain assembly-only. Fusion support counts are not
assembly EM abundance estimates.

### Python API

```python
from fin.pipeline import PipelineConfig, PipelineRunner

config = PipelineConfig.from_profile(
    "real-drna",
    bam_path="reads.bam",
    gtf_path="annotations.gtf",
    genome_fasta_path="genome.fa",
    fastq_path="reads.fq",
    signal_path="signal.blow5",
    signal_format="slow5",
    work_dir="./pyfin_output",
    output_gtf="./pyfin_output/assembly.gtf",
    output_tsv="./pyfin_output/scores.tsv",
    fusion_enabled=False,
    use_gpu=True,
    em_sigma=1.0,
)
config.validate()  # resolves profile routes and effective compatibility state

runner = PipelineRunner(config)
runner.setup()
results = runner.run()
runner.cleanup()

for cid, qr in sorted(results.items(), key=lambda x: -x[1].abundance):
    print(f"{qr.candidate_id}\t{qr.abundance:.2f}\t{qr.num_assigned_reads}\t{qr.source}")

from fin.fusion.pipeline import run_fusion_analysis

fusion_results = run_fusion_analysis(
    "reads.bam", "genome.fa", "./pyfin_output",
    gtf_path="annotations.gtf", sequence_path="reads.fq",
)
```

## Per-Transcript Scoring

Signal scoring is layered into read-to-transcript assignment: M1 supplies the
mappy structural mask and M2 refines exact ties with per-event junction NLL.
These feed the assignment engine, whose soft responsibilities yield each
candidate's `abundance`, `num_reads`, and `tpm` in `scores.tsv`.

The `coherence_score`, `discrimination_score`, and `combined_score` columns in
`scores.tsv` are retained for output-schema stability but are always `0.0`: they
were produced by a separate composite scorer that has been removed. Do not rely
on them.

## Fusion Detection

Enabled with `--fusion`, after ordinary assembly finishes. The independent path:

1. Scans BAM records and SA tags without relying on assembly intervals. It
   normalizes forward-query coordinates, including reverse and hard-clipped
   alignments, and checks optional FASTQ/FASTA sequences for conflicts.
2. Realigns full molecules and clipped arms to the genome, retaining normal
   full-length and paralogous competitors rather than choosing a single partner.
3. Preserves ordered N-arm paths. A-B-C requires reads spanning the whole chain;
   independent A-B and B-C reads cannot manufacture it. The default minimum is
   two distinct read IDs; supplementary records do not add votes.
4. Proposes internal splice and fusion-junction corrections within
   `--fusion-junction-tolerance` (default 3 bp). Read-window sequence evidence,
   both-side anchors and complete stitched-transcript validation are required.
   Canonical motifs and annotation are nominations, not evidence of expression.
5. Clusters coordinate tuples with bounded diameter `--max-dist` (default 6 bp).
   Microhomology/near-tied coordinates remain explicit alternatives, even when
   the ordered gene chain is supported. Readthrough candidates are classified
   separately; missing annotation produces loci, not fabricated gene names.

After assembly, `--fusion` reads the final `assembly.gtf` and
`junction_support.tsv`. Final-model junctions with at least three indel-clean
reads supply nearby nominations and structural annotations, never extra fusion
read votes. Fusion-edge donor/acceptor support is not support for the edge itself.
These aggregated counts are not claimed to be independently signal-validated.

Fusion also uses the supplied signal with assembly's krill eventalign and local
**summed-LLR** reducer, pore, GPU request, window flank and decision margin.
Only already sequence-supported complete paths are eligible: conflicting
coordinate alternatives are removed, while identical-template ambiguity,
missing signal, incomplete contrasts and weak coverage explicitly abstain.
If all alternatives conflict, the read no longer supplies a fusion vote.
Abstention retains clearly labeled sequence-only evidence, not a signal pass.
This shares the scoring scheme, not assembly's AS-tie gate, candidates or EM.
The standalone API remains sequence-only unless `signal_path` is provided.
These RNA-level calls do not establish genomic rearrangements.

`fusion/fusions.tsv` contains supported, unresolved-coordinate and insufficient-
support chains with explicit statuses; do not count every row as an accepted
call. `fusions.bedpe` has N-1 rows per chain with chain ID and edge index, using
representative coordinates. `fusion_transcripts.fa` and `fusion_arms.jsonl`
retain alternative templates and exon paths; no fictitious multi-chromosome GTF
is emitted. `fusion_reads.jsonl`, `junction_corrections.jsonl`, `junction_signal.jsonl` and
`fusion_manifest.json` preserve decisions, limits, settings and provenance.
The TSV's `signal_status` distinguishes full/partial local signal support from
sequence-only calls. See `docs/FUSION_SIGNAL_VALIDATION.md` for scope and checks.

**Fused-transcript level.** On top of the chains, `fusion_models.tsv` reports
one row per *fused transcript model*: which annotated isoform of each gene
took part (`compatible_transcript_ids`, per-arm `arm_match_classes`), the
conventional junction name from exon ranks (`gA:e13>gB:e2`), the breakpoints,
and molecule support split into `support_unique` (fits this model alone),
`support_compatible` (fits this and other models -- a 5'-degraded molecule
usually cannot single out one isoform) and `support_partial`. Multi-gene
models (A-B-C) appear here with three arms exactly when complete molecules
span the whole chain. Two FASTAs are written and never share a record id:
`fusion_models.fa` is evidence-constrained (the envelope of what molecules
observed, reference-projected), `fusion_models.annotated.fa` extends the
outer arms to compatible annotated ends and says so in the header
(`kind=annotated_completion extension_bases_5p=... extension_observed=false`).
`fusion_model_exons.tsv` is the multi-chromosome exon table (0-based
half-open genomic and transcript offsets, per-exon provenance) that
reconstructs either FASTA; `fusion_model_evidence.jsonl` carries the read ids;
`fusion_model_collapsed.jsonl` lists observed structures that were folded into
a longer model as possible 5' truncations, with the models they could be. See
`docs/FUSION_TRANSCRIPT_MODELS_2026-09-15.md`.

Publication switches the `fusion/` symlink atomically between run-local
directories. Check `fusion_manifest.json` for `status="complete"` before using
results. A failed rerun exposes a failed manifest, not the previous calls.

Migration from the older prototype: the public `--fusion` path no longer
mixes fusion candidates into assembly EM, writes root-level `fusions.bedpe`, or
uses 500-bp clustering/flank defaults. The legacy `PipelineConfig.fusion_*`
fields remain for internal compatibility; use `FusionSettings` and
`run_fusion_analysis` for the isolated Python API. Existing un-managed
`fusion/` directories must be moved aside before publication.

The independent scan covers the input BAM, not just assembly's scheduled
intervals or read-limit subset. Use an explicitly subsetted BAM when a restricted
fusion analysis is intended; candidate/search caps abstain or fail rather than
silently sample a prefix. Implementation acceptance and the synthetic/real-signal
assembly byte comparisons are recorded in `docs/FUSION_VALIDATION_2026-09-12.md`.

### Using individual modules

```python
# Candidate discovery only
from fin.candidates import discover_candidates
from fin.io.interval_manager import GenomicInterval
from fin.io.io_gtf import GTFReader

interval = GenomicInterval(chrom="chr1", start=1000, end=50000, strand="+")
with GTFReader("annotations.gtf") as gtf:
    gtf.parse()
    genome_seq = open("genome.fa").read()  # simplified
    candidates = discover_candidates(interval, "reads.bam", gtf, genome_seq)

# Eventalign parsing only
from fin.scoring import parse_eventalign_tsv, build_distance_matrix
scores = parse_eventalign_tsv("eventalign.tsv", candidate_lengths={"tx1": 1000})
dist = build_distance_matrix(scores, read_ids=["r1", "r2"], candidate_ids=["tx1", "tx2"])

# EM assignment only
from fin.analysis.assignments import em_with_coherence
R, assignments, lls = em_with_coherence(
    dist_read_to_tx, dist_read_to_read,
    sigma=1.0, beta=0.5, use_gpu=True
)
```

## GPU Acceleration

PyFIN uses GPU acceleration in three areas:

| Component | Technology | Fallback |
|-----------|-----------|----------|
| DTW pairwise distances | CUDA C++ extension (`fin._dtw`) | scipy/numpy CPU |
| EM algorithm (matrix ops) | CuPy (`cupy.matmul`, `cupy.exp`) | numpy CPU |
| DTW variable-length pairs | CUDA `dtw_distance()` per pair | numpy row-vectorized |

### GPU memory safety

The DTW pairwise batch function checks available GPU memory before allocation. For large datasets (>10,000 reads), it will raise a `MemoryError` with a descriptive message rather than crashing with a CUDA OOM error.

### Disabling GPU

```python
# Pipeline level
config = PipelineConfig(..., use_gpu=False)

# EM algorithm level
R, assignments, lls = em_with_coherence(
    dist_read_to_tx, dist_read_to_read,
    use_gpu=False  # Force numpy CPU
)
```

## Module Structure

```
fin/
  candidates/               # Transcript candidate discovery
    dataclasses.py          # IntronChain, TranscriptCandidate, CandidateSet
    intron_chains.py        # CIGAR -> intron chain, 3' clustering
    discovery.py            # Full candidate discovery per interval
  scoring/                   # External tool scoring
    external_tools.py       # minimap2/f5c subprocess wrappers
    eventalign_parser.py    # f5c TSV -> distance matrix
    signal_dtw.py           # Signal extraction + pairwise DTW
  pipeline/                  # Orchestration
    config.py               # PipelineConfig dataclass
    runner.py               # PipelineRunner (Phase 1-5)
  analysis/                  # Statistical methods
    assignments.py          # EM with coherence (CuPy GPU)
    quantification.py       # Probability-weighted abundance
    clustering.py           # 3' end position clustering
  io/                        # I/O for bioinformatics formats
    io_bam.py               # BAM reader (pysam)
    io_gtf.py               # GTF reader
    io_pod5.py              # POD5 signal reader
    io_slow5.py             # SLOW5/BLOW5 signal reader
    io_fasta.py             # FASTA reader
    io_fastq.py             # FASTQ reader
    interval_manager.py     # Genomic interval management
  _dtw/                      # CUDA DTW extension
    __init__.py             # Python API (dtw_distance, dtw_pairwise)
    dtw_api.cpp             # C++/CUDA bindings
    dtw.hpp                 # CUDA DTW kernel
  utils/                     # Utilities
    log_config.py           # Logging configuration
    sequences.py            # Sequence utilities
```

## Configuration Reference

### PipelineConfig parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `bam_path` | str | required | Path to BAM file |
| `gtf_path` | str | None | Path to GTF annotation file |
| `genome_fasta_path` | str | required | Path to genome FASTA |
| `fastq_path` | str | required | Path to reads FASTQ |
| `signal_path` | str | required | Path to signal file (SLOW5/BLOW5/POD5) |
| `signal_format` | str | "slow5" | Signal format: "slow5" or "pod5" |
| `work_dir` | str | "./pyfin_work" | Working directory for intermediate files |
| `output_gtf` | str | None | Output assembly GTF path |
| `output_tsv` | str | None | Output per-candidate scoring TSV path |
| `output_bedpe` | str | None | Legacy coupled-fusion output; not used by public `--fusion` |
| `three_prime_threshold` | int | 24 | 3' end clustering distance (bp) |
| `isoform_fraction_locus` | str | "family" | Use persistent splice families for minor-isoform denominators; "overlap" restores legacy behavior |
| `post_selection_refit` | bool | False | Renormalize beta=0 M2 responsibilities over final survivors; named profiles enable it |
| `complete_chain_rescue_reads` | int | 0 | P2 evidence rescue (`--complete-chain-rescue-reads`): a multi-exon candidate whose weakest junction is cleanly attested by enough exact complete-chain reads (gtf >= N; novel multi-junction >= N+1; novel single-junction >= 10) is exempt from the m2 support gate and soft-mass filter. SIRV4 expressed-truth evidence at N=2 on top of the strict-positive gate: F1 89.02 -> 89.29, `c_trunc3p50` +1.96, 3/156 cells regress <= 0.42. Opt-in; 0 disables |
| `junction_election_mode` | str | "off" | K-way signal-ballot junction election (`--junction-election-mode off\|audit\|require`): nearby junction variants with a cigar-level support imbalance are put to a per-read raw-signal vote (krill eventalign NLL over the contested windows; quorum 5, supermajority 0.75, window purity, event/margin safeguards). `audit` writes junction_election.tsv and is byte-inert. Decided elections pick the correct junction 97.6% of the time on SIRV truth. Requires m2_em + `--signal` + families discovery. Opt-in |
| `junction_election_burden` | str | "prosecution" | Burden of proof for `require`: `prosecution` retains unless the vote convicts (conservative; SIRV F1 89.29 -> 90.38, 0 cells down); `defense` drops nominated junctions unless the vote exonerates (aggressive; SIRV F1 89.29 -> **91.87**, top 87 -> **112/156**, worst cell -0.06, 73 junctions exonerated). Defense is SIRV-fit; re-validate on real dRNA before use |
| `junction_election_ballot` | str | "signal" | Development ablation: `signal`/`m1`/`screen`/`shuffle` ballots (only `signal` is a production setting) |
| `assignment_channels` | str | "none" (profiles `real-drna*`: "ends") | Endpoint evidence in the m2_em assignment posterior (`--assignment-channels none\|tes\|ends`). `tes`: the sample's fitted 3'/polyA kernel scores every polyA-PASS read against each tied model's claimed TES; `ends` adds the one-sided 5' start term (a read cannot begin upstream of the model it came from). An annotated model whose ends the reads contradict (a polyA read past its TES, or -- under `ends` -- a read starting upstream of its TSS, in neither case explained by another annotated model on the same or a wobble chain) gets a read-extent sibling `<id>.ext` with the same intron chain; selection keeps one of the pair by posterior mass. `tes` applies this to single-exon annotation matched by overlap, `ends` also to multi-exon annotated keys. Kernels must be calibrated (>= 50 ends, narrow scale <= 25 bp) or the channel abstains; without `--signal` the 3'/polyA term abstains (no calls = no evidence) and only the 5' term applies, and that term needs a `--gtf` (it scores read starts against the annotated TSS), so an unguided run without signal has no endpoint evidence at all. Requires m2_em. Chosen by cross-dataset maximin: +1.54 / +0.64 / +1.59 F1 (SIRV4 / Sequin / HEYA8, expression contract) with Sn and Pr both up on every grid |
| `folded_tss_recovery_mode` | str | "off" (profiles `real-drna*`: "require") | Post-EM recovery of generation-folded 5'-truncation shadows: each folded shadow's TSS is tested against the parent's read-5'-end hazard profile; `audit` writes tss_shadow_recovery.tsv only, `require` re-emits a supported shadow as its own candidate and reruns the interval's assignment once. Maximin: +0.05 / 0 / +0.03 F1, never negative |
| `containment_evidence_mode` | str | "off" | Containment decision layer (`--containment-evidence-mode off\|audit\|require`): contained short-in-long candidate pairs are routed through degradation-aware TES-cluster / termination-hazard tests (tri-state verdicts; abstention always retains). `audit` writes containment_evidence.tsv and is byte-inert; `require` drops contradicted shorts (needs an independently supported parent) and exempts supported ones from the competition gates. Requires m2_em. Opt-in |
| `require_gtf_positive_mass` | bool | False | Strict-positive GTF evidence gate (`--require-gtf-positive-mass`): a GTF-annotated transcript must end selection with EM soft-mass > 0 or it is dropped as unexpressed, even when every abundance floor is 0 (enforced both at the pre-refit A4 floor and again after the survivor refit). Never adds catalog entries; novel/fusion semantics unchanged. Requires `gtf_path`. Opt-in; no profile enables it. SIRV4 expressed-truth evidence (with `--min-gtf-abundance 0 --no-floor-gtf-abundance`): F1 85.58 -> 89.02, zero regressing cells |
| `em_sigma` | float | 1.0 | EM temperature (lower = harder assignments) |
| `em_max_iter` | int | 1000 | Maximum EM iterations |
| `em_tol` | float | 1e-4 | EM convergence tolerance |
| `use_gpu` | bool | True | Enable GPU acceleration (DTW + EM) |
| `max_reads_per_interval_for_dtw` | int | 2000 | DTW subsampling cap per interval |
| `fusion_enabled` | bool | False | Legacy interval-coupled API; keep False and use `run_fusion_analysis` for isolation |
| `fusion_min_support` | int | 2 | Legacy coupled-fusion support; independent API uses `FusionSettings.min_support` |
| `fusion_max_dist` | int | 500 | Legacy value; independent API defaults to a 6-bp bounded cluster diameter |
| `fusion_flank_bp` | int | 500 | Legacy value; independent API defaults to 35-bp read-window flanks |
| `f5c_path` | str | "f5c" | Path to f5c binary |
| `max_reads` | int | None | Limit total reads processed |

## Testing

### Run all unit tests

```bash
# With the correct Python environment (needs pysam, numpy, scipy)
pytest tests/unit/ -v

# Run specific test modules
pytest tests/unit/test_intron_chains.py -v      # Intron chain extraction
pytest tests/unit/test_eventalign_parser.py -v   # Eventalign TSV parsing
pytest tests/unit/test_quantification.py -v      # Quantification
pytest tests/unit/test_assignments.py -v         # EM algorithm
pytest tests/unit/test_dtw.py -v                 # DTW (CUDA tests skipped without GPU)
```

### Test with GPU

```bash
# DTW tests require CUDA GPU
pytest tests/unit/test_dtw.py -v  # CUDA tests auto-skip if no GPU

# Verify GPU is being used
python -c "
from fin._dtw import is_available
from fin.analysis.assignments import CUPY_AVAILABLE
print(f'CUDA DTW: {is_available()}')
print(f'CuPy EM:  {CUPY_AVAILABLE}')
"
```

### Test data

Test data is located in `tests/testdata/`:
- `RNA004.test.pod5` - POD5 signal file
- `test.fa` - Reference FASTA
- `RNA004.test.bam` - Aligned reads

## Development

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run tests with coverage
pytest tests/unit/ --cov=fin --cov-report=term-missing

# Lint
ruff check fin/
```

## Benchmarking

The published comparison is produced by the runners under
`experiments/prod_validation/`, not by a generic harness:

| Step | Script |
|---|---|
| Run one (sample, annotation-condition) cell of the grid | `experiments/prod_validation/sirv4_grid_runner.py` |
| Score assemblies against the expression truth | `experiments/prod_validation/expression_rescore.py` |
| Build the per-cell tool-vs-tool table | `experiments/prod_validation/gw_headtohead.py` |

The truth contract is the **expression** contract: a transcript counts as
truth for a sample when it is in the reference annotation AND NanoCount
reports TPM > 0 for that same sample. Structural agreement is gffcompare
transcript level `=` (exact intron chain). Sensitivity and precision are
always reported together.

`benchmarks/` holds analysis and plotting utilities used while developing
those runs. It is not a one-command reproduction of the paper.

> **Reproduction recipe: not yet written.** These runners read dataset trees
> configured for the authors' machine (`sirv4_grid_runner.py`), so a reader
> cannot currently re-run the evaluation from a clean checkout. A release-time
> recipe -- input manifest with checksums, annotation-condition generation,
> competitor commands and versions, and the table-generation step -- is
> outstanding. This is stated rather than implied: previously this section
> advertised a `run_benchmark.sh` that wrote `"status": "ran"` markers without
> ever invoking the tools. That script has been removed.

## License

MIT. Full text in [LICENSE](LICENSE); Copyright (c) 2026 Yuelong Chen.

Third-party components keep their own terms. In particular the signal backend
`krill` is installed separately from its own index and is not vendored here.
