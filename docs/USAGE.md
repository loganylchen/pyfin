# Using PyFIN

Installation is in the root [`README.md`](../README.md). **Read the krill
section there first** -- the shipped profile requires it and will not start
without it.

This document covers the three things installation does not: preparing inputs
that actually match each other, diagnosing a run that finished without doing
what you expected, and reading the numbers that come out.

For what the pipeline does internally, see
[`CURRENT_ALGORITHM.md`](CURRENT_ALGORITHM.md).

---

## 1. Preparing inputs

PyFIN takes five inputs. The README names them; this section is about getting
them to agree with each other, which is where most first runs fail.

| Input | Flag | Notes |
|---|---|---|
| Alignments | `--bam` | coordinate-sorted, indexed, **spliced** |
| Genome | `--genome` | FASTA, same build and contig names as the BAM |
| Annotation | `--gtf` | optional; same build and contig names again |
| Reads | `--fastq` | the same reads the BAM contains |
| Raw signal | `--signal` | SLOW5/BLOW5 or POD5, same reads again |

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
`chr1` with a GTF on `1` will run and silently annotate nothing, because no
interval will ever match a gene.

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

### 1.3 Reads and signal must be the same molecules

Three identifiers have to line up:

- the BAM's `query_name`,
- the FASTQ header -- PyFIN takes **the first whitespace-separated token**
  after `@`, so `@abc-123 runid=... ch=42` is read id `abc-123`,
- the read id inside the raw-signal file.

A quick overlap check before a long run:

```bash
samtools view sample.sorted.bam | cut -f1 | sort -u > /tmp/bam.ids
zcat reads.fq.gz | awk 'NR%4==1{print substr($1,2)}' | sort -u > /tmp/fq.ids
comm -12 /tmp/bam.ids /tmp/fq.ids | wc -l     # compare against UNIQUE bam ids:
wc -l < /tmp/bam.ids
```

Compare against the number of **unique query names**, not the alignment-record
count -- a read with supplementary alignments contributes several records.

The signal file needs the same check and the command depends on your format
(`slow5tools skim`, or a POD5 reader). Doing it is worth the minutes: a
signal file from a different run is the most common cause of a run that
completes with signal contributing nothing.

If any overlap is near zero you probably basecalled and aligned different
files, or the headers carry a prefix.

### 1.4 Signal format and pore model are not auto-detected

```bash
--signal reads.blow5                 # default: --signal-format slow5
--signal reads.pod5 --signal-format pod5
```

`--signal-format` defaults to `slow5`. Passing a POD5 without the flag is a
common first-run failure.

The pore model defaults to `rna002` (`fin/pipeline/config.py:749`) and is
**not** detected from the file. A mismatched model degrades signal scoring
rather than refusing to run.

**There is currently no CLI flag for it.** The CLI does not populate
`krill_pore` when it builds `PipelineConfig`, so changing it means using the
Python API and setting the field on the config object directly. If your data
is RNA004 chemistry, this is a real limitation to plan around rather than a
flag to pass.

---

## 2. Diagnosing a run

### 2.1 "It finished, but did signal do anything?"

The import check in the README tells you krill is *installed*. It does not
tell you signal scoring *contributed*. These are different, and a run can
complete with signal abstaining on every read.

Signal can legitimately abstain, per read, when:

- there were no exact alignment-score ties to refine -- nothing for the
  junction term to break;
- no valid scored contrast could be obtained at all;
- the backend initialised but produced no usable aligner, or a batch failed
  and the per-read fallback also produced nothing.

What is **not** abstention: a structure selected on a sufficient score margin
*without* a verified covered span. Coverage controls whether that structure
also contributes to the prior (`fin/pipeline/assignment.py:1127`); its absence
is not a failure to score.

None of these stop the run, and none of them mean your inputs are wrong. What
*is* worth investigating is abstention on essentially everything, which
usually means the signal file and the BAM are not the same reads (see §1.3).

The closest direct evidence is in the interval log: the per-interval `ties`
and `refined` counters (`fin/pipeline/runner.py:2654`). `ties` is how many
exact alignment-score ties were found; `refined` counts ties where **at least
two hypotheses received valid scores** -- note that this increments *before*
the structure-selection margin is applied
(`fin/pipeline/assignment.py:597`, `:1127`), so it means "a junction contrast
was scored", not "the assignment changed".

Large `ties` with `refined` at zero everywhere is worth investigating, but it
does **not** by itself identify mismatched ids. Valid, matching inputs can
also produce no disagreement window, invalid candidate payloads, too few
scored hypotheses, or a backend failure. Neither counter is a quality metric.

The manifest tells you the configuration, not whether signal contributed:

```bash
# what actually ran, with the resolved configuration
python -m json.tool pyfin_out/run_manifest.json | head -40
```

`run_manifest.json` records the **resolved** profile and every override. When
a document and a run disagree about a setting, the manifest is right.

### 2.2 "It failed with ModuleNotFoundError: krill"

Expected, and it means exactly what it says: the shipped profile resolves
`m2_metric = "summed_llr"`, and the assignment layer imports krill on that
path before checking whether any read needs refinement. Install krill (root
README) or run with `--m2-metric off` if you deliberately want a
sequence-only result.

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
(`fin/pipeline/assignment.py:390`, `fin/scoring/krill_aligner.py:66`). So the
two accelerations are separate installs but not separate switches.

The CuPy check in the README is only about the assignment/EM step.

---

## 3. Reading the outputs

Default assembly outputs land in `--output-dir`:

| File | What it is |
|---|---|
| `assembly.gtf` | the final transcript set, with `abundance` as an attribute |
| `scores.tsv` | one row per reported candidate: abundance, read counts, TPM, diagnostics |
| `abundance_refit.json` | mass accounting for the refit -- orphaned mass, conservation |
| `run_manifest.json` | the resolved configuration, source hash, environment |
| `junction_support.tsv` | observed junction evidence per interval |

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
