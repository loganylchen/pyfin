# Changelog

What changed for someone who runs PyFIN. Mechanism is in
[`docs/CURRENT_ALGORITHM.md`](docs/CURRENT_ALGORITHM.md); benchmark numbers are
not recorded here.

## Unreleased (to become 1.0.0)

### A reported value changed

- **The `tpm` column of `scores.tsv` now defaults to the direct-RNA
  convention**, `1e6 * abundance / sum(abundances)`. It used to divide by
  transcript length first, the short-read convention. Pass
  `--tpm-length-normalize` to get the old column. The assembly
  (`assembly.gtf`, which models are reported, their abundance) is identical
  either way; only this column -- in `scores.tsv` and, when it is written, in
  `scores.unfiltered.tsv` -- and the TPM figures inside `abundance_refit.json`
  differ.

No other change in this release is meant to alter what is assembled. What was
checked is under "Verification" below; what it does not cover is said there.

### Fixed

- **An extreme-depth interval no longer exhausts memory.** The signal
  tie-break handed every tied read of an interval to the signal backend in one
  call and kept all of its output. The mitochondrial genome of a deep library
  is one interval of more than 10^5 reads, and a run on it was killed at
  480 GB. The reads are now handed over in chunks that close at
  `m2_batch_pairs` (2000) read-candidate pairs. If a call to the backend fails,
  only the reads of that chunk fall back to one-by-one scoring; before, the
  whole interval did.
- **The `pyfin` command is installed.** Only `fin` was; every document used
  `pyfin`.
- **`pyslow5 >= 1.0` is required.** Older versions were allowed by the package
  metadata and silently dropped every read's signal.
- Messages logged by the command itself (the resolved profile, the run
  summary) were not shown when PyFIN was started as `python -m fin.cli`.

### Added

- **Input checks before the run starts** (`--preflight`, on by default). A
  BAM without an index, a genome on different contig names or of a different
  build, an annotation that shares no contig name with the BAM, a signal file
  that holds none of the BAM's reads, a missing signal backend: each stops the
  run at once with the fix, instead of failing hours in or finishing with a
  result that only looks like biology. Partial mismatches (some contigs
  missing from the genome, some sampled reads without signal) are warnings.
  The output directory is not created until the checks pass. The checks cannot
  tell an annotation of another assembly build whose contig names match.
  `--no-preflight` skips them.
- **`--pore {auto,rna002,rna004}`.** RNA004 could not be selected from the
  command line before. `auto` (the default) reads the sequencing kit, or
  failing that the sampling rate, from the signal file. RNA004 support is
  **experimental**: see [`docs/RNA004_STATUS.md`](docs/RNA004_STATUS.md).
- **`--signal-format auto`** (the default) takes the format from the file name.
- **`run_summary.json`**: how many reads reached assignment, how many were in
  an alignment-score tie, how many of those ties the signal scored, how many
  models came out, the time and the peak memory of one process. It answers
  "did the signal score any junction?" without reading a log. It is written
  last, and an older one is removed when a run starts writing, so a run that
  failed part-way leaves none.
- `run_manifest.json` records the PyFIN version and the versions of krill,
  mappy, pysam, pyslow5, pod5, numpy and scipy.
- **`--version`**, printing the same versions.
- **`--help-all`.** `--help` now lists about thirty options, the ones a run
  needs; the experimental switches and diagnostics (more than a hundred) are
  under `--help-all`. They work exactly as before.
- A small real dataset with raw signal, `tests/data/sirv1_mini`, and
  end-to-end tests that run the whole command on it.

### Changed

- **`--fastq` is optional.** The assembly takes read sequences from the BAM
  and never read the FASTQ. With `--fusion` it is still used, to recover the
  sequence of hard-clipped supplementary alignments.
- **POD5 input is refused, with the conversion command.** The signal backend
  reads SLOW5/BLOW5 only; handed a POD5 it crashed. Convert once with
  `blue-crab p2s`.
- The output directory no longer receives one empty sub-directory per
  interval (a whole-genome run left about twelve thousand).
- **No more `logs/` directory.** Importing PyFIN, and every worker process of
  a run, used to create `./logs/fin_<timestamp>.log` in whatever directory it
  was started from; a whole-genome campaign accumulated gigabytes of them, and
  PyFIN could not start in a directory it could not write to. PyFIN now logs
  to the console only. Redirect it to keep a log
  (`pyfin ... > run.log 2>&1`).
- Python 3.10 or newer, on Linux x86-64: the signal backend ships wheels for
  CPython 3.10 to 3.12 only.
- The package is pure Python. `setup.py` no longer builds a CUDA extension;
  nothing had used it since DTW moved to the signal backend.

### Known limitations

These come from our whole-genome benchmark on fifteen SG-NEx RNA002 libraries
(truth: reference transcripts expressed in each library; the report is not
part of this repository yet). They are observations on that benchmark, not
general rankings.

- **No annotation at all (`de novo`).** Against StringTie3, recall was equal
  and precision about four points lower, in every library.
- **Deep libraries.** The precision difference to StringTie3 fell as the
  number of reads in a library rose, under every annotation condition tested.
- **RNA004** runs, and is not yet benchmarked.
- **Speed.** A Python loop calls the aligner (mappy) once for every read and
  every candidate of its interval. A whole-genome library takes hours on eight
  workers.
- **Memory with `--threads N`** is roughly N times the memory of one worker.
- **A CUDA build of the signal backend** chooses its CPU or GPU kernel per
  batch, so with it the chunking above, like the composition of an interval,
  can move a score. The shipped images use the CPU backend.

### Verification

- The test suite, including the end-to-end tests on real signal.
- A byte-level baseline of a small assembly (`assembly.gtf`, `scores.tsv`,
  `abundance_refit.json`) recorded before this series: with
  `--tpm-length-normalize` all three files still match it exactly.
- For the chunked signal scoring: on three loci of a real run, the previous
  code and the new code give the same `assembly.gtf` and `scores.tsv`, at the
  default chunk size and with chunks forced down to 1 and to 7 pairs.

These are observations with one set of library versions (krill 0.1.5) and
with every call to the signal backend succeeding. They do not cover a run in
which a backend call fails and reads fall back to one-by-one scoring, nor a
CUDA build of the backend.
