# Using PyFIN

Installation is in the root [`README.md`](../README.md). **Read the krill
section there first** -- the shipped profile requires it and will not start
without it.

This document covers the three things installation does not: preparing inputs
that actually match each other, diagnosing a run that finished without doing
what you expected, and reading the numbers that come out.

PyFIN checks most of what section 1 describes by itself, before it starts
(section 1.5). The section still matters: the checks tell you *that* two
inputs disagree, and this is where to read *why* they must agree.

For what the pipeline does internally, see
[`CURRENT_ALGORITHM.md`](CURRENT_ALGORITHM.md).

---

## 1. Preparing inputs

PyFIN takes three inputs and two optional ones. The README names them; this
section is about getting them to agree with each other, which is where most
first runs fail.

| Input | Flag | Notes |
|---|---|---|
| Alignments | `--bam` | coordinate-sorted, indexed, **spliced**; primary records must carry the read sequence |
| Genome | `--genome` | FASTA, same build and contig names as the BAM |
| Raw signal | `--signal` | SLOW5/BLOW5, the same reads again |
| Annotation | `--gtf` | optional; same build and contig names again |
| Reads | `--fastq` | optional; the assembly takes read sequences from the BAM. Only `--fusion` reads it |

### 1.1 The BAM must be spliced, sorted and indexed

Align direct-RNA reads with a splice-aware preset. `minimap2 -ax splice -uf
-k14` is the usual choice for direct RNA; `-uf` forces the forward transcript
strand, which is what direct RNA gives you.

```bash
minimap2 -ax splice -uf -k14 -t 8 genome.fa reads.fq.gz \
  | samtools sort -@ 4 -o sample.sorted.bam -
samtools index sample.sorted.bam
```

The index is not optional: the reader opens `<bam>.bai` and fails if it is
absent. A missing index and a genuinely unreadable BAM produce different
errors -- if `samtools quickcheck sample.sorted.bam` passes and PyFIN still
cannot open it, check for the `.bai` before anything else.

An unspliced alignment will not error. It will produce reads with no introns,
so the pipeline finds no multi-exon structure and you get a mono-exon-heavy
result that looks like a biology problem and is not.

### 1.2 The genome, the BAM and the GTF must be the same build

Not just the same assembly version -- the same **contig naming**. A BAM on
`chr1` with a GTF on `1` can never match a gene. PyFIN refuses to start when
the two share no contig name at all; with `--no-preflight`, or through the
Python API, such a run goes ahead and silently annotates nothing.

```bash
# a naming SPOT-CHECK -- not proof the assemblies agree
samtools idxstats sample.sorted.bam | cut -f1 | head -3
grep -v '^#' annotation.gtf | cut -f1 | sort -u | head -3
grep '^>' genome.fa | cut -c2- | cut -d' ' -f1 | head -3
```

Comparing the first few names tells you whether the `chr` prefix convention
matches. It does not establish that the three are the same assembly build, nor
that every contig in one exists in the others.

If they disagree, fix the annotation or the reference -- do not rename
contigs in the BAM header alone, because the genome FASTA is read by contig
name too.

### 1.3 The alignment and the signal must be the same molecules

Two identifiers have to line up: the BAM's `query_name` and the read id inside
the raw-signal file. (With `--fusion --fastq`, the FASTQ header is a third:
PyFIN takes **the first whitespace-separated token** after `@`, so
`@abc-123 runid=... ch=42` is read id `abc-123`.)

PyFIN checks this before it starts: it looks up the first 200 reads of the BAM
in the signal file. If none is there it searches the whole BAM for one that
is, and refuses to start only when the two files share no read at all -- the
usual causes being a signal file from another run, or read ids that were
rewritten (a prefix, a suffix).

A signal file that holds only part of the BAM's reads is accepted: reads
without signal are assembled from sequence alone. PyFIN warns when some of the
200 sampled reads have no signal. It does not measure coverage beyond that
sample -- if the first reads all have signal and later ones do not, nothing
is reported at start-up; `run_summary.json` shows afterwards how many ties
the signal scored.

### 1.4 Signal format and pore model

```bash
--signal reads.blow5                 # the format is taken from the file name
--signal reads.dat --signal-format slow5
```

**POD5 is not read.** PyFIN's signal backend opens the signal file itself and
reads SLOW5/BLOW5 only. Convert a POD5 run once, losslessly:

```bash
blue-crab p2s pod5_dir/ -o reads.blow5       # pip install blue-crab
```

Given a POD5, PyFIN stops with that command rather than start.

**The pore model** is chosen by `--pore {auto,rna002,rna004}`. `auto`, the
default, reads every read group's header in the signal file: the sequencing
kit, or failing that the sampling rate (3012 Hz for RNA002, 4000 Hz for
RNA004). What it does when the header is not clear-cut:

| The signal file says | `--pore auto` |
|---|---|
| one chemistry | uses it |
| nothing (no kit, no known sampling rate) | uses rna002, with a warning |
| a chemistry in some read groups, nothing in others | uses that chemistry, with a warning |
| different chemistries in different read groups | refuses: pass `--pore` yourself |
| a DNA experiment in every read group | refuses |

A wrong pore model does not fail. It scores the signal against the wrong
current levels and quietly degrades the result, which is why PyFIN never
guesses between two. RNA004 is selectable and not yet benchmarked; see
[`RNA004_STATUS.md`](RNA004_STATUS.md).

### 1.5 What is checked before a run starts

Unless `--no-preflight` is given, PyFIN verifies the following before it
creates the output directory. An error stops the run; a warning does not.
"Sampled" means the first 200 mapped primary records of the BAM. The checks
write nothing, with one exception: a genome FASTA without a `.fai` index gets
one, next to the FASTA.

| Check | Outcome when it fails |
|---|---|
| the signal backend imports (unless `--m2-metric off`) | error, with the install command |
| the signal file is SLOW5/BLOW5 | error, with the conversion command (also under `--no-preflight`) |
| the BAM opens and has an index | error: `samtools index` |
| the sampled primary records carry a read sequence | error if none of them does |
| some sampled read is spliced | warning if none is, among at least 50 sampled (section 1.1); says so if the reference names look like transcripts |
| contigs with mapped reads exist in the genome | error if none does (naming); warning if some do not |
| those contigs have the same length in the BAM header and the genome | error: different builds |
| the annotation shares a contig name with the BAM | error (section 1.2) |
| the annotation has `exon` records | warning |
| the signal file holds reads of the BAM | error if it holds none of them; warning if some of the sampled reads are missing (section 1.3) |
| the pore model is unambiguous | see section 1.4 |
| the signal backend can estimate poly(A) tails for that pore | warning: the run completes without poly(A) calls (krill before 0.1.6, on RNA004). Every read abstains on the 3' endpoint evidence; where enabled, endpoint refinement uses end sharpness alone and the polyA+5' filter (off in the default profile) is skipped |

What it cannot check: that a GTF whose contig names match is the *same
assembly build* as the alignment (section 2.3); how much of the BAM beyond the
sample has a sequence, a spliced alignment or signal; or that your reads are
what you think they are.

---

## 2. Diagnosing a run

### 2.1 "It finished, but did signal do anything?"

Look at `run_summary.json`:

```json
"signal": {
  "scoring": "summed_llr",
  "pore": "rna002",
  "tied_reads": 17,
  "scored_contrasts": 5
}
```

- `tied_reads` -- reads whose best alignment score is shared by two or more
  candidates. These are the only reads the signal is asked about; a read with
  a single best candidate is already decided.
- `scored_contrasts` -- ties in which **at least two candidates received a
  valid signal score**. It is counted (`fin/pipeline/assignment.py:644`,
  `:684`) before the decision margin is applied (`:1181`), so it means "a
  junction contrast was scored", not "the assignment changed".

These two counters are about the **junction** scoring only. The poly(A)
evidence that the endpoint channels use is also read from the raw signal, and
is not counted here.

The same two counters are in the log, per interval, as `ties` and `refined`.

Signal can legitimately abstain, per read, when:

- there were no exact alignment-score ties to refine -- nothing for the
  junction term to break;
- no valid scored contrast could be obtained at all;
- the backend initialised but produced no usable aligner, or a batch failed
  and the per-read fallback also produced nothing.

What is **not** abstention: a structure selected on a sufficient score margin
*without* a verified covered span. Coverage controls whether that structure
also contributes to the prior (`fin/pipeline/assignment.py:1180`); its absence
is not a failure to score.

None of these stop the run, and none of them mean your inputs are wrong. What
*is* worth investigating is many `tied_reads` with `scored_contrasts` at zero.
It does **not** by itself identify mismatched ids -- the start-up check
(section 1.5) has already established that the signal file holds reads of the
BAM. Valid, matching inputs can also produce no disagreement window, invalid
candidate payloads, too few scored hypotheses, or a backend failure. Neither
counter is a quality metric.

`run_summary.json` also records the models written, the elapsed time, the
peak memory of one process and any warning the start-up checks raised. It is
written last, and a summary already in the directory is removed when a run
starts writing -- so a run that failed part-way leaves no summary. (A run that
the start-up checks refused never started: whatever an earlier run left in the
directory, its summary included, is untouched.) With a quantification engine
other than the default (`--quant-mode`), the two signal counters are `null`,
which means "not recorded", not zero.

The manifest tells you the configuration, not whether signal contributed:

```bash
# what actually ran, with the resolved configuration
python -m json.tool pyfin_out/run_manifest.json | head -40
```

`run_manifest.json` records the **resolved** profile and every override, and
the versions of PyFIN, krill, mappy, pysam, pyslow5, pod5, numpy and scipy.
When a document and a run disagree about a setting, the manifest is right.

### 2.2 "It stopped: the krill backend cannot be imported"

It means exactly what it says: the shipped profile resolves
`m2_metric = "summed_llr"`, which scores the raw signal with krill. PyFIN
checks for it at start-up and prints the install command (root README). Note
that the package of the same name on PyPI is a different project.
`--m2-metric off` turns the junction scoring off, which lets a run start
without krill. That run is sequence-only; with krill installed, `--m2-metric
off` still leaves the poly(A) evidence, which is read from the signal too.

Through the Python API, or with `--no-preflight`, the same condition surfaces
later, as `ModuleNotFoundError` when the first interval reaches assignment.

### 2.3 "Nothing was annotated"

Check §1.2 -- but note what the check can and cannot tell you. Annotation
lookup matches on contig **name** and coordinate overlap, not on assembly
identity (`fin/io/io_gtf.py:493`), and matching records become `source="gtf"`
candidates (`fin/candidates/discovery.py:563`).

So a contig-name mismatch can prevent annotation matching entirely, while a
**wrong-build GTF whose contig names happen to match will still supply
annotated candidates** -- at coordinates that belong to a different assembly.
An all-novel result is not guaranteed by a mismatch, and a completed run does
not establish that your inputs agree.

### 2.4 GPU questions

`--no-gpu` forces the CPU path for the EM **and is passed to krill as well**.
A CUDA-built krill requests the GPU when `krill.built_with_cuda` is true and a
device is visible, and falls back to CPU if initialisation fails
(`fin/pipeline/assignment.py:410`, `fin/scoring/krill_aligner.py:66`). So the
two accelerations are separate installs but not separate switches.

CuPy concerns the assignment/EM step only; krill needs no CuPy.

---

## 3. Reading the outputs

Default assembly outputs land in `--output-dir`:

| File | What it is |
|---|---|
| `assembly.gtf` | the final transcript set, with `abundance` as an attribute |
| `scores.tsv` | one row per reported candidate: abundance, read counts, TPM, diagnostics |
| `run_summary.json` | what the run did: models, reads in a tie, signal contrasts scored, time, peak memory, warnings (section 2.1) |
| `run_manifest.json` | what the run was asked to do: the resolved configuration, PyFIN and library versions, source hash, environment |
| `abundance_refit.json` | mass accounting for the refit -- orphaned mass, conservation |
| `junction_support.tsv` | observed junction evidence per interval |
| `tss_shadow_recovery.tsv` | the folded 5'-truncated chains tested for a start of their own; written only when there were any |

`--fusion` adds a `fusion/` directory with its own tables and FASTAs; see
[`README.md`](../README.md) and the fusion documents in
[`README.md` of this directory](README.md).

### 3.1 abundance, num_reads and confidence are three different things

- **`abundance`** is fractional. It is the soft, EM-assigned mass: a read
  compatible with three isoforms contributes a fraction to each. This is the
  quantity to use.
- **`num_reads`** is the hard count -- reads whose best assignment was this
  candidate, ties broken deterministically (so a tied read still lands
  somewhere). It is an integer, and it is **not** simply "abundance rounded
  down": ten reads splitting 0.6/0.4 give soft masses 6 and 4 but hard counts
  10 and 0. A candidate's hard count can exceed its soft mass. Under the
  refit, hard counts cover the reads assigned to survivors, not the whole
  input library.
- **`confidence`** summarises how concentrated a candidate's assigned
  responsibility was. **It is not a probability that the isoform exists.** A
  candidate with one unambiguous read can have high confidence; that says the
  read was unambiguous, not that the transcript is real.

### 3.2 TPM carries a convention you should choose deliberately

By default the `tpm` column is the **direct-RNA convention**:
`tpm = 1e6 * abundance / sum(abundances)`. Direct RNA is not fragmented: one
molecule yields at most one read regardless of its length, so the read count
already is the per-molecule quantity. This is what NanoCount and
`salmon --ont` report. Note that it is a million-sum share, not the raw
abundance -- the two coincide only when the abundances happen to sum to a
million.

`--tpm-length-normalize` divides abundance by transcript length first. That is
the short-read convention, which exists because fragmentation makes a longer
transcript yield proportionally more reads. It was PyFIN's default before
v1.0; pass the flag to reproduce older tables.

Neither is an unbiased estimate of molarity. Against the known input
concentrations of the Sequin spike-in the default agrees better overall, and
the length-normalised value agrees better on transcripts of 2 kb and longer.
Pick the one that matches what you are comparing against, and say which you
used.

**This flag changes reporting only.** No structural selection decision reads
the TPM column; the transcript set is identical either way.

### 3.3 The diagnostics in `abundance_refit.json`

Three numbers are worth looking at:

- **alignment-unassigned reads** -- reads that reached assignment but matched
  no candidate well enough. A large share suggests the candidate set is
  missing real structure, or the reads are from something not in the
  reference.
- **orphaned vs rescued mass** -- these are different and the refit reports
  both (`fin/analysis/abundance_refit.py:371`). A read that still has a
  surviving candidate is renormalised over the survivors: that is **rescued**
  mass. A read whose candidates were *all* dropped has nowhere to go and
  stays **unassigned** -- orphaned mass is not redistributed, it is accounted
  for (`fin/analysis/abundance_refit.py:238`, `:251`). A large orphaned figure
  means selection removed everything some reads matched.
- **conservation error** -- should be at or near zero. A non-trivial value
  means mass was lost or duplicated in the refit, which is a bug worth
  reporting, not a property of your data.

None of the three is a quality gate. They are there so you can see what the
pipeline did with your reads.
