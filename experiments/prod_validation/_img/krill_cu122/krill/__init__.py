"""krill - in-memory direct-RNA eventalign (RNA002 / RNA004), CPU + GPU.

Given a read's own basecalled sequence and its raw signal, return eventalign-format
results as a list of per-read dicts of numpy arrays.

Example
-------
    import krill
    aligner = krill.Aligner(pore="rna004", use_gpu=False)
    result = aligner.align({
        "read_id": "read1",
        "sequence": "ACGUACGU...",          # A/C/G/T/U, any case (U->T); other bases -> ValueError
        "signal": raw_signal,                # np.ndarray, raw ADC samples (float32)
        "digitisation": 8192.0,
        "offset": 5.0,
        "range": 1437.0,
        "sample_rate": 4000.0,
    })
    # result is a list[dict]; one dict per read with per-column numpy arrays.
"""

import numbers

from ._krill import Aligner, built_with_cuda, gpu_device_count

# DTW pairwise-distance primitive (additive downstream of the aligner; see dtw.py)
from .dtw import (
    dtw_distance,
    dtw_pairwise,
    dtw_pairwise_varlen,
    dtw_multi_position_pairwise,
    backend as dtw_backend,
    is_available as dtw_available,
    extension_available as dtw_extension_available,
    cleanup as dtw_cleanup,
)


def _is_int(x):
    """True for Python ints and numpy integers (e.g. np.int64), but not bool."""
    return isinstance(x, numbers.Integral) and not isinstance(x, bool)

__all__ = ["Aligner", "align_blow5", "align_read_variants",
           "align_reads_variants",
           "built_with_cuda", "gpu_device_count",
           "dtw_distance", "dtw_pairwise", "dtw_pairwise_varlen",
           "dtw_multi_position_pairwise",
           "dtw_backend", "dtw_available", "dtw_extension_available",
           "dtw_cleanup"]


def align_blow5(blow5_path, fastq_path, pore="rna004", use_gpu=None,
                batch_size=0, num_thread=0, cuda_dev_id=0):
    """Optional helper: read reads from a fastq + a slow5/blow5 file and align them.

    Requires `pyslow5` and a fastq reader (`pyfastx`). This is a convenience path
    on top of the pure in-memory API; the core engine itself does not need files.
    """
    import pyslow5
    import pyfastx

    aligner = Aligner(pore=pore, use_gpu=use_gpu, batch_size=batch_size,
                      num_thread=num_thread, cuda_dev_id=cuda_dev_id)

    s5 = pyslow5.Open(blow5_path, "r")
    try:
        reads = []
        for name, seq, _ in pyfastx.Fastq(fastq_path, build_index=False):
            rec = s5.get_read(name, pA=False, aux=["read_id"])
            if rec is None:  # read not present in the blow5; silently skipped
                continue
            reads.append({
                "read_id": name,
                "sequence": seq,
                "signal": rec["signal"].astype("float32"),
                "digitisation": float(rec["digitisation"]),
                "offset": float(rec["offset"]),
                "range": float(rec["range"]),
                "sample_rate": float(rec["sampling_rate"]),
            })
    finally:
        s5.close()

    return aligner.align(reads)


def align_read_variants(blow5_path, read_id, sequences, pore="rna004",
                        use_gpu=None, num_thread=0, cuda_dev_id=0,
                        aligner=None, start=None):
    """Align ONE read's signal against MANY candidate sequences.

    Fetches the raw signal for `read_id` from a slow5/blow5 file (pyslow5 builds
    the .idx automatically on first access if it is missing), then runs
    eventalign for each candidate sequence against that same signal.

    Parameters
    ----------
    blow5_path : str
        Path to a .slow5 / .blow5 file containing `read_id`.
    read_id : str
        Name of the read whose signal is used.
    sequences : str | dict | iterable
        The candidate sequence(s) to align against the signal. Accepts:
          * a single sequence string, or
          * a list/tuple of sequence strings, or
          * a dict {label: sequence} to tag each result with a label.
    pore, use_gpu, num_thread, cuda_dev_id :
        Forwarded to `Aligner` (ignored if `aligner` is supplied).
    aligner : Aligner, optional
        Reuse an existing Aligner instead of constructing one (useful when
        calling this repeatedly to avoid re-loading the pore model / GPU).
    start : int | dict | iterable, optional
        Coordinate offset added to the `position` column of each result, for
        when a candidate sequence is a *truncated* slice of a transcript and
        you want positions in transcript coordinates. The offset is the
        transcript coordinate of the candidate's first base, so the first
        emitted position becomes `start` instead of 0. Accepts:
          * a single int applied to every candidate, or
          * a dict {label: offset} keyed by the same labels as `sequences`
            (missing labels default to 0), or
          * a list/tuple of ints aligned with `sequences` by position.
        Defaults to 0 for all candidates.

    Returns
    -------
    list[dict]
        One eventalign result dict per candidate sequence, in input order.
        Each dict has the usual per-column numpy arrays plus a "variant_label"
        key identifying which candidate it came from.

    Example
    -------
        import krill
        res = krill.align_read_variants(
            "reads.blow5", "read-uuid",
            {"ENS1": "ACGUACGU...", "ENS2": "ACGCACGU..."},
            start={"ENS1": 1, "ENS2": 5},
            pore="rna004", use_gpu=True)
        for r in res:
            print(r["variant_label"], r["position"][0])   # -> ENS1 1 / ENS2 5
    """
    import pyslow5

    # normalise `sequences` -> list of (label, seq)
    if isinstance(sequences, str):
        variants = [("0", sequences)]
    elif isinstance(sequences, dict):
        variants = [(str(k), v) for k, v in sequences.items()]
    else:
        variants = [(str(i), s) for i, s in enumerate(sequences)]
    if not variants:
        raise ValueError("`sequences` is empty")

    # normalise `start` -> per-variant offset aligned with `variants`
    if start is None:
        offsets = [0] * len(variants)
    elif isinstance(start, dict):
        offsets = [int(start.get(label, 0)) for label, _ in variants]
    elif _is_int(start):
        offsets = [int(start)] * len(variants)
    else:
        offsets = [int(s) for s in start]
        if len(offsets) != len(variants):
            raise ValueError("`start` length (%d) != number of sequences (%d)"
                             % (len(offsets), len(variants)))

    # pyslow5 auto-creates <blow5_path>.idx on first get_read if absent
    s5 = pyslow5.Open(blow5_path, "r")
    try:
        rec = s5.get_read(read_id, pA=False)
    finally:
        s5.close()
    if rec is None:
        raise KeyError("read_id %r not found in %s" % (read_id, blow5_path))

    signal = rec["signal"].astype("float32")
    common = dict(
        signal=signal,
        digitisation=float(rec["digitisation"]),
        offset=float(rec["offset"]),
        range=float(rec["range"]),
        sample_rate=float(rec["sampling_rate"]),
    )

    reads = [dict(read_id=read_id, sequence=seq, start=off, **common)
             for (_, seq), off in zip(variants, offsets)]

    own = aligner is None
    if own:
        aligner = Aligner(pore=pore, use_gpu=use_gpu,
                          batch_size=max(len(reads), 1),
                          num_thread=num_thread, cuda_dev_id=cuda_dev_id)

    results = aligner.align(reads)
    for (label, _), r in zip(variants, results):
        r["variant_label"] = label
    return results


def align_reads_variants(blow5_path, reads_variants, pore="rna004",
                         use_gpu=None, num_thread=0, cuda_dev_id=0,
                         aligner=None, start=None, batch_size=0):
    """Align MANY reads, each against MANY candidate sequences, in one batch.

    The multi-read counterpart of `align_read_variants`. Instead of one Python
    call per read (each a tiny batch), this flattens every (read_id, candidate)
    pair into a SINGLE `aligner.align(...)` call so the C++ engine's `num_thread`
    workers parallelise across the whole workload — removing the per-call Python
    / GIL / result-copy overhead that dominates when looping over thousands of
    reads with a handful of candidates each.

    Parameters
    ----------
    blow5_path : str
        Path to a .slow5 / .blow5 file containing every read_id (pyslow5 builds
        the .idx automatically on first access if missing).
    reads_variants : dict | iterable
        The reads and their candidate sequences. Accepts:
          * a dict {read_id: sequences}, or
          * an iterable of (read_id, sequences) pairs,
        where each `sequences` is itself anything `align_read_variants` accepts:
        a single string, a list/tuple of strings, or a {label: sequence} dict.
    pore, use_gpu, num_thread, cuda_dev_id :
        Forwarded to `Aligner` (ignored if `aligner` is supplied).
    aligner : Aligner, optional
        Reuse an existing Aligner instead of constructing one.
    start : dict, optional
        Per-read coordinate offsets: {read_id: start_spec}, where each
        `start_spec` matches the `start=` argument of `align_read_variants`
        (a single int, a {label: offset} dict, or a list aligned with that
        read's candidates). Reads not present default to offset 0.
    batch_size : int, optional
        If >0, submit the flattened pairs to `align()` in chunks of at most this
        many pairs (bounds peak Python-side memory for very large workloads).
        The engine still re-chunks internally by the Aligner's own batch_size;
        this only limits how many pairs are materialised per `align()` call.
        Default 0 = submit everything in one call.

    Returns
    -------
    dict
        {read_id: list[dict]} — for each read, the list of per-candidate result
        dicts in that read's candidate order, each carrying a "variant_label"
        key (and "read_id"). Reads whose signal is missing from the file are
        omitted.

    Example
    -------
        import krill
        out = krill.align_reads_variants(
            "reads.blow5",
            {"read-A": {"ref": "ACGU...", "alt": "ACGC..."},
             "read-B": ["GGAC...", "GGAU..."]},
            pore="rna004", num_thread=8)
        for rid, res in out.items():
            for r in res:
                print(rid, r["variant_label"], r["status"])
    """
    import pyslow5

    # normalise reads_variants -> list of (read_id, [(label, seq), ...])
    if isinstance(reads_variants, dict):
        items = list(reads_variants.items())
    else:
        items = list(reads_variants)

    def norm_seqs(sequences):
        if isinstance(sequences, str):
            return [("0", sequences)]
        if isinstance(sequences, dict):
            return [(str(k), v) for k, v in sequences.items()]
        return [(str(i), s) for i, s in enumerate(sequences)]

    def norm_starts(start_spec, variants):
        if start_spec is None:
            return [0] * len(variants)
        if isinstance(start_spec, dict):
            return [int(start_spec.get(label, 0)) for label, _ in variants]
        if _is_int(start_spec):
            return [int(start_spec)] * len(variants)
        offs = [int(s) for s in start_spec]
        if len(offs) != len(variants):
            raise ValueError("`start` length (%d) != number of sequences (%d)"
                             % (len(offs), len(variants)))
        return offs

    s5 = pyslow5.Open(blow5_path, "r")

    reads = []          # flattened align() input
    pair_meta = []      # parallel: (read_id, label) for each entry in `reads`
    order = []          # read_ids in input order (for stable output ordering)
    seen = set()

    try:
        for read_id, sequences in items:
            variants = norm_seqs(sequences)
            if not variants:
                continue
            rec = s5.get_read(read_id, pA=False)
            if rec is None:  # signal missing from the file; this read is dropped
                continue
            signal = rec["signal"].astype("float32")
            common = dict(
                signal=signal,
                digitisation=float(rec["digitisation"]),
                offset=float(rec["offset"]),
                range=float(rec["range"]),
                sample_rate=float(rec["sampling_rate"]),
            )
            offsets = norm_starts(start.get(read_id) if isinstance(start, dict)
                                  else start, variants)
            for (label, seq), off in zip(variants, offsets):
                reads.append(dict(read_id=read_id, sequence=seq, label=label,
                                  start=off, **common))
                pair_meta.append((read_id, label))
            if read_id not in seen:
                seen.add(read_id)
                order.append(read_id)
    finally:
        s5.close()

    if not reads:
        return {}

    own = aligner is None
    if own:
        aligner = Aligner(pore=pore, use_gpu=use_gpu,
                          batch_size=max(len(reads), 1),
                          num_thread=num_thread, cuda_dev_id=cuda_dev_id)

    # run, optionally chunked to bound peak Python-side memory
    results = []
    if batch_size and batch_size > 0:
        for base in range(0, len(reads), batch_size):
            chunk = aligner.align(reads[base:base + batch_size])
            # the engine numbers read_index from 0 within each align() call; add
            # the chunk base back so read_index stays a global input position.
            for r in chunk:
                r["read_index"] = base + int(r["read_index"])
            results.extend(chunk)
    else:
        results = aligner.align(reads)

    # regroup by read_id, preserving per-read candidate order
    out = {rid: [] for rid in order}
    for (rid, label), r in zip(pair_meta, results):
        r["variant_label"] = label
        r["read_id"] = rid
        out[rid].append(r)
    return out
