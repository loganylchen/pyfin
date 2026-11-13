# RNA004: what is in place, and what is not

RNA004 can be selected and runs through the same code as RNA002. **It has not
been benchmarked.** Treat an RNA004 assembly as experimental until it has.

This page separates what has been checked from what is assumed. It contains no
accuracy numbers because none exist yet.

## How to run

```bash
pyfin --bam sample.bam --genome genome.fa --signal sample.blow5 \
      --gtf annotation.gtf --output-dir pyfin_out/
```

RNA004 needs krill 0.1.6 or later (see "Poly(A) estimation" below). That is
the version the images and the install instructions pin; an older
installation is upgraded with:

```bash
pip install --upgrade krill --no-deps \
    --index-url https://loganylchen.github.io/krill-dist/simple/
```

Nothing RNA004-specific is needed on the command line. `--pore auto`, the
default, reads the sequencing kit from the signal file (`sqk-rna004`), or
failing that the sampling rate (4000 Hz), and selects the RNA004 pore model.
The start-up line says what it chose and why:

```
preflight ok: signal=slow5, pore=rna004 (sequencing kit 'sqk-rna004'), ...
```

and `run_summary.json` records it under `signal.pore`. Pass `--pore rna004`
to choose it yourself, for instance for a converted file that lost its header
fields; the cases in which `auto` refuses or warns are in
[`USAGE.md`](USAGE.md) section 1.4.

A POD5 run has to be converted to BLOW5 first (`blue-crab p2s`), as for
RNA002.

## What has been checked

- **The pore model reaches every signal call.** One configuration field,
  `krill_pore`, is passed to the junction scoring, the poly(A) estimation,
  and the fusion analysis's signal check (`fin/pipeline/runner.py`,
  `fin/pipeline/assignment.py`, `fin/scoring/polya.py`, `fin/cli.py`). Before
  v1.0 it existed but the command line never set it, so RNA004 could only be
  reached through the Python API.
- **Detection.** On a real RNA004 BLOW5 (kit `sqk-rna004`, 4000 Hz) `auto`
  resolves to `rna004`; on all fifteen SG-NEx RNA002 libraries it resolves to
  `rna002`. A wrong pore model does not fail -- it degrades signal scoring --
  so `auto` never picks between two chemistries found in one file.
- **The signal backend aligns it.** krill 0.1.5 accepts `pore="rna004"`
  (a 9-mer model; RNA002 is a 5-mer model) for eventalign, which is what the
  junction scoring uses.
- **Poly(A) estimation needs krill 0.1.6.** krill 0.1.5 has no RNA004
  poly(A) model and refuses to build the estimator; PyFIN then runs with no
  poly(A) calls at all. Every read abstains on the 3' endpoint evidence (the
  5' evidence, which comes from the alignments, is unaffected); where they
  are enabled, endpoint refinement judges 3' ends by end sharpness alone, and
  the polyA+5' filter -- off in the default profile -- is skipped rather than
  applied. The start-up checks say so, with the upgrade command. krill 0.1.6
  added the model: on forty reads of an RNA004 HeLa library it returned an
  estimate for thirty-nine (thirty-six passing its own quality flag), where
  0.1.5 returned none. Whether those lengths are accurate was not checked
  here.
- **krill 0.1.6 against 0.1.5 on RNA002.** The two versions give the same
  assembly on two whole-genome benchmark runs (one H9 library, with the full
  annotation and with none) and the same `assembly.gtf` and `scores.tsv` on
  three loci of another run. The benchmark was run with 0.1.5 and the other
  benchmark runs have not been repeated; the Docker images, CI and the
  install instructions use 0.1.6.
- **One end-to-end run on real RNA004 data.** 3,939 reads of a HeLa native
  RNA004 library (4,013 sequenced), aligned to GRCh38, with GENCODE v44 and
  `--pore auto`. The run completes in about five minutes on four workers
  with a peak of 3.4 GB, the pore is read from the kit, and a second run
  reproduces the assembly. With krill 0.1.6 it reports 888 models (882
  annotated, 6 novel); with 0.1.5, which has no poly(A) calls, 880, a
  different set. That the two differ is expected. Which is closer to the
  truth is not known.
- **The signal reaches the scoring windows on RNA004.** In that run 1,663
  reads were in an alignment-score tie. For all but 61 of them the signal
  was not consulted, by design: the tied candidates had different numbers of
  introns (1,190 attempts) or no junction on which they differ (416
  attempts); four reads have two attempt records. For those 61, the signal
  backend aligned 137 of 151 read-candidate pairs, and every event of those
  137 had a usable score and mapped to a genomic position. Seven pairs had
  a scoring window overlapping the genomic span of their aligned events, and
  all seven had events in a window; the other 130 had neither. The same
  pattern holds on an RNA002 locus (43 pairs with an overlapping window, all
  with events in one; 8 with neither). Three reads ended with a scored
  contrast and one of
  them was decided (margin 4.06, over 25 events per structure); the other
  two tied exactly. That is far too few to say anything about how well the
  signal discriminates on RNA004.
- **The junction-window arithmetic assumes no k-mer size.** The junction
  score is computed from events inside a window of +/- 8 genomic bp around
  each splice boundary on which the candidates differ. krill reports every
  event at the centre base of its k-mer, for either model, so the window
  arithmetic does not depend on the k-mer length. One 5-mer assumption does
  remain in the code: an event the backend could not model is recognised by
  the model k-mer `NNNNN`. According to the backend's source such events come
  only from its HMM re-decoding, which the default junction scoring does not
  request; the match-fraction tie-break (`fin/scoring/krill_tiebreak.py`),
  which is off by default and has no command-line option, does request it and
  has not been examined for RNA004.

## What has not been established

- **Accuracy.** No RNA004 library has been assembled against a truth set.
  The run above shows that the pipeline works on RNA004 input, not that its
  output is right.
- **The junction decision threshold.** The signal is consulted for a read
  whose best alignment score is shared by several candidates with the same
  number of introns, when at least two of their structures have scored events
  in the +/- 8 bp windows; otherwise it abstains. Such a read is assigned to
  one structure when the summed log-likelihood ratio between the two best
  structures is at least 1.0. Both numbers were chosen on RNA002.
  RNA004 current levels have a different noise scale and a 9-mer context, so
  the same threshold may be more or less permissive there. This is the first
  thing a benchmark has to look at.
- **The endpoint evidence.** The 3'/poly(A) and 5' kernels are fitted on each
  sample, so they adapt to the data in principle; that they behave on RNA004
  as on RNA002 has not been examined. With krill 0.1.5 the poly(A) part is
  absent altogether (above).
- **A test on real RNA004 signal.** The bundled example and the end-to-end
  tests that run are RNA002. The start-up behaviour for RNA004 (detecting the
  pore, the warning above) is tested on synthetic headers and a stand-in
  backend. Two older integration tests (`test_cpu_fallback.py`,
  `test_discovery_gtf_vs_no_gtf.py`) are written against an RNA004 dataset,
  `tests/testdata/RNA004.test.*`, that is not in the repository; they are
  skipped.

## What a benchmark will need

The same protocol as the RNA002 benchmark: reference transcripts expressed in
each library as truth, exact intron-chain matches, recall and precision
reported together, the annotation supplied at several levels of completeness
and corruption, the same competing assemblers. The two RNA002-derived
constants above must not be re-tuned on the libraries they are then judged
on; if they need recalibrating, that takes data held out from the comparison.
