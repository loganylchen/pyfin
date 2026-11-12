# PyFIN - signal-aware isoform assembly for nanopore direct RNA

PyFIN assembles full-length transcript models from nanopore **direct RNA**
reads aligned to a genome. It discovers known and novel isoforms from the
reads' own intron chains, uses the **raw current signal** to assign reads that
sequence alignment alone cannot place, treats a reference annotation as
evidence to be confirmed rather than as truth, and reports an abundance for
every model it keeps. An optional, independent analysis reports multi-gene
fusion transcripts.

## Where to read what

| You want to | Go to |
|---|---|
| Install it | [Installation](#installation) below |
| Check the installation and see a run | [Quick start](#quick-start-the-bundled-example) below |
| Run your own data | [Running your own data](#running-your-own-data) below |
| Prepare inputs, diagnose a run, interpret the outputs | [`docs/USAGE.md`](docs/USAGE.md) |
| Know what the algorithm actually does | [`docs/CURRENT_ALGORITHM.md`](docs/CURRENT_ALGORITHM.md) |
| Know what changed, and what is known not to work well | [`CHANGELOG.md`](CHANGELOG.md) |
| Use RNA004 data | [`docs/RNA004_STATUS.md`](docs/RNA004_STATUS.md) |
| Find any other document | [`docs/README.md`](docs/README.md) -- the index, including which documents are historical and what is known to be missing |

Several documents in this repository read as authoritative because each was,
on the day it was written. `docs/README.md` says which one is current for
which question; `ALGORITHM.md` and `PIPELINE_AS_IMPLEMENTED.md` in particular
are dated snapshots and carry banners saying so.

## Pipeline overview

```
genome BAM + genome FASTA + raw signal (SLOW5/BLOW5)   [+ GTF, optional]
                |
        Input checks (before the output directory is created)
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
                  - run_summary.json   (what the run did)
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

**Requirements.** Linux on x86-64, Python 3.10 to 3.12, and a C compiler with
the zlib headers (`mappy`, the minimap2 binding, is built from source; on
Debian or Ubuntu: `apt-get install build-essential zlib1g-dev`). The supported
platform is set by the signal backend, which ships wheels for exactly that.

### 1. PyFIN

```bash
git clone https://github.com/loganylchen/pyfin.git
cd pyfin
pip install .
```

This installs the `pyfin` command (`fin` is the same command). `mappy` and
`pysam` come in as dependencies -- PyFIN calls them in-process and never
shells out to the `minimap2` or `samtools` programs.

### 2. krill, the signal backend

```bash
pip install "krill==0.1.5" --no-deps \
    --index-url https://loganylchen.github.io/krill-dist/simple/
```

**Use that index.** An unrelated package called `krill` exists on PyPI;
`pip install krill` without the index installs the wrong thing. For the same
reason krill is not, and cannot be, a declared dependency of PyFIN.

krill 0.1.5 is the version every published benchmark number was produced with.
It supplies the eventalign that scores the raw signal. The shipped profile
needs it: without it PyFIN stops at start-up and prints this install command.
`--m2-metric off` turns the junction signal scoring off and lets a run start
without krill; an installed krill is then still used for poly(A) evidence.

For RNA004 data install krill 0.1.6 or later instead: 0.1.5 has no RNA004
poly(A) model ([`docs/RNA004_STATUS.md`](docs/RNA004_STATUS.md)).

### 3. Optional: GPU for the assignment step

```bash
pip install ".[gpu]"      # CuPy
```

CuPy accelerates the read-to-transcript assignment matrices. It is separate
from krill, and both obey `--gpu/--no-gpu`. The Docker images use the CPU
build of krill on purpose: its eventalign calls are many and small, and
across several worker processes the CPU build is faster.

### Docker

```bash
docker build -f Dockerfile.cpu -t pyfin:cpu .      # or Dockerfile.gpu
docker run --rm -v "$PWD":/data pyfin:cpu pyfin --version
```

### Check the installation

```bash
pyfin --version            # PyFIN and the libraries that decide a result
scripts/smoke_test.sh      # runs the bundled example; a few seconds
```

## Quick start: the bundled example

`tests/data/sirv1_mini` is one locus of a real direct-RNA run: 91 reads, their
raw signal, the reference sequence and eight annotated transcripts.

```bash
pyfin \
    --bam    tests/data/sirv1_mini/reads.bam \
    --genome tests/data/sirv1_mini/reference.fa \
    --signal tests/data/sirv1_mini/reads.blow5 \
    --gtf    tests/data/sirv1_mini/annotation.gtf \
    --output-dir sirv1_out --no-gpu
```

It ends with a line like this one, which is also what `run_summary.json` holds:

```
Run summary: 10 models (8 annotated, 2 novel) from 2 of 2 intervals; 17 reads in an
alignment-score tie, 5 signal contrasts scored (summed_llr, rna002); ...
```

"17 reads in an alignment-score tie" are reads that two or more candidate
transcripts explain equally well by sequence; "5 signal contrasts scored" are
the ties in which the raw signal could be compared against the candidates. A
run in which that number is zero scored no splice junction with the signal.
(The exact counts depend on the krill version; these are from 0.1.5.)

## Running your own data

| Input | Option | Notes |
|---|---|---|
| Alignments | `--bam` | reads aligned to the **genome** with a splice-aware preset, coordinate-sorted and indexed |
| Genome | `--genome` | the FASTA the reads were aligned to |
| Raw signal | `--signal` | SLOW5 or BLOW5 of the same reads |
| Annotation | `--gtf` | optional; same assembly and contig names as the BAM |
| Reads | `--fastq` | optional; only `--fusion` uses it |

```bash
# align (direct RNA: forward strand, splice-aware)
minimap2 -ax splice -uf -k14 -t 8 genome.fa reads.fq.gz \
  | samtools sort -@ 4 -o sample.bam - && samtools index sample.bam

# POD5 must be converted once; PyFIN's signal backend reads SLOW5/BLOW5
blue-crab p2s pod5_dir/ -o sample.blow5        # pip install blue-crab

pyfin --bam sample.bam --genome genome.fa --signal sample.blow5 \
      --gtf annotation.gtf --output-dir pyfin_out/ --threads 8
```

**Before it starts, PyFIN checks that the inputs agree with each other**: the
BAM index, contig names and lengths against the genome, a contig name shared
with the annotation, that the signal file holds reads of the BAM, and that the
signal backend is installed. A mismatch it can prove stops the run at once and
says how to fix it; a partial one (some contigs missing from the genome, some
reads without signal) is a warning. The output directory is not created until
the checks pass. What they cannot establish is listed in
[`docs/USAGE.md`](docs/USAGE.md) section 1.5 -- in particular, an annotation
of another assembly build whose contig names happen to match is not detected.
`--no-preflight` skips the checks.

**Pore model.** `--pore auto` (the default) reads the sequencing kit from the
signal file: RNA002 and RNA004 are recognised. RNA004 runs and is not yet
benchmarked -- see [`docs/RNA004_STATUS.md`](docs/RNA004_STATUS.md).

**Workers and memory.** `--threads N` runs N worker processes, one genomic
interval each. Memory grows with N: a deep human library needs on the order of
100-150 GB with eight workers. `run_summary.json` reports the peak memory of
one process.

**Options.** `pyfin --help` lists the options an ordinary run needs.
`pyfin --help-all` lists every option. The additional ones are of two kinds:
settings the profile owns (several are on in the shipped profile;
[`docs/CURRENT_ALGORITHM.md`](docs/CURRENT_ALGORITHM.md) section 0 lists
them), and experimental switches and diagnostics that no profile enables. An
ordinary run needs none of them.

### Operating profiles

The CLI resolves a named scientific operating point before constructing
`PipelineConfig`. The product default is balanced `--profile real-drna`; use
`--profile real-drna-precision` for a higher read-support floor on novel models
(it has not been benchmarked in combination with the current defaults),
`--profile sirv` for synthetic SIRV benchmarking, and `--profile custom` for raw
option defaults. Any explicitly supplied CLI option overrides its profile value.

| Profile | Intended data | Key resolved behavior |
| --- | --- | --- |
| `real-drna` | Biological nanopore dRNA, balanced | Tight summed LLR; abundance floors novel > 1, annotated >= 1; novel mono >=5 hard reads; read-supported junction consensus; final survivor-set abundance refit |
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

### Outputs

- `assembly.gtf` -- the assembled transcripts, each with its abundance,
  confidence, hard-read count and source (`gtf` or `novel`)
- `scores.tsv` -- one row per model: abundance, read count and `tpm`
  (by default the direct-RNA convention, a share of abundance per million;
  see [`docs/USAGE.md`](docs/USAGE.md) section 3.2)
- `run_summary.json` -- what the run did: models, reads in a tie, signal
  contrasts scored, elapsed time, peak memory of one process, warnings
- `run_manifest.json` -- what the run was asked to do: resolved profile and
  configuration, PyFIN and library versions, source hash
- `junction_support.tsv` -- read support of every junction observed
- `abundance_refit.json` -- mass conservation of the final abundance re-fit
- `tss_shadow_recovery.tsv` -- the folded 5'-truncated chains that were tested
  for a start of their own (written only when there were any)
- with `--fusion`: `fusion/` (see [Fusion detection](#fusion-detection))

Optional diagnostics (`--candidate-evidence`, `--endpoint-refine`,
`--catalog-mode`, `--m2-contrast-stats-jsonl`, `--selection-audit`,
`--write-unfiltered-scores`) write further files; each is described by
`pyfin --help-all`.

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

### Assembly plus fusion detection

```bash
pyfin \
    --bam reads.bam --gtf annotations.gtf --genome genome.fa \
    --signal signal.blow5 --fastq reads.fq \
    --output-dir pyfin_out/ \
    --fusion --min-support 2 --max-dist 6 --flank-bp 35
```

Adds results under `pyfin_out/fusion/`. Ordinary `assembly.gtf`, `scores.tsv`
and `abundance_refit.json` remain assembly-only. Fusion support counts are not
assembly EM abundance estimates. `--fastq` lets the fusion analysis recover the
sequence of hard-clipped supplementary alignments.

### Python API

The command line is the supported interface. The same run through the API:

```python
from fin.pipeline import PipelineConfig, PipelineRunner
from fin.preflight import run_preflight

checks = run_preflight(                      # raises PreflightError on a mismatch
    bam_path="reads.bam", genome_path="genome.fa", signal_path="signal.blow5",
    signal_format="auto", m2_metric="summed_llr", gtf_path="annotations.gtf",
)
config = PipelineConfig.from_profile(
    "real-drna",
    bam_path="reads.bam",
    genome_fasta_path="genome.fa",
    signal_path="signal.blow5",
    gtf_path="annotations.gtf",
    signal_format=checks.signal_format,
    krill_pore=checks.pore,                  # rna002 or rna004, read from the signal
    work_dir="./pyfin_output",
    output_gtf="./pyfin_output/assembly.gtf",
    output_tsv="./pyfin_output/scores.tsv",
    use_gpu=False,
)
config.validate()  # resolves profile routes and effective compatibility state

runner = PipelineRunner(config)
runner.setup()
try:
    models = runner.run()
finally:
    runner.cleanup()

for model in sorted(models.values(), key=lambda m: -m.abundance):
    print(f"{model.candidate_id}\t{model.abundance:.2f}\t{model.num_assigned_reads}\t{model.source}")
```

Importing `fin` configures no logging and writes no file; call
`fin.setup_logger("fin", level="INFO")` to see the run's messages.

```python
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

## Testing and development

```bash
pip install -e ".[dev]"        # pytest, ruff, build
pip install "krill==0.1.5" --no-deps \
    --index-url https://loganylchen.github.io/krill-dist/simple/

pytest                         # tests/unit and tests/integration
ruff check --select E9,F63,F7,F82 fin tests/unit tests/integration
```

The suite includes end-to-end tests that run the whole command on the bundled
signal fixture (`tests/integration/test_sirv1_end_to_end.py`); they are skipped
when krill is not installed, and CI refuses to run without it.

Continuous integration (`.github/workflows/`) runs the suite on Python 3.10 to
3.12, builds the wheel and runs the smoke test from a clean environment, and,
for releases, builds the images and runs the smoke test inside them before
publishing.

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

## Known limitations

They are listed in [`CHANGELOG.md`](CHANGELOG.md#known-limitations): assembly
without any annotation, very deep libraries, RNA004, speed, and memory with
many workers.

## License

MIT. Full text in [LICENSE](LICENSE); Copyright (c) 2026 Yuelong Chen.

Third-party components keep their own terms. In particular the signal backend
`krill` is installed separately from its own index and is not vendored here.
