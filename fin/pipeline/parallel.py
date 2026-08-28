"""Interval-level process parallelism for the assembly pipeline.

The per-interval work (:meth:`PipelineRunner.process_interval`) is wholly
self-contained and the cross-interval aggregation
(:func:`aggregate_across_intervals`) is commutative, so intervals can be run in
independent worker processes and merged in any order.

Concurrency model (see ``fin/pipeline/runner.py`` ``run()``):
  * ``spawn`` start method -- the parent has already imported numpy, so a forked
    worker would inherit the parent's BLAS thread pool and the per-worker thread
    caps set in ``cli.main()`` would not take effect; spawn re-imports fresh.
    CUDA / pysam / Pod5 handles are also not fork-safe.
  * Static GPU/CPU partition -- exactly ``gpu_workers`` of the ``threads`` workers
    build a GPU-backed runner (``config.use_gpu = True``); the rest are CPU-only.
    This bounds live CUDA contexts to ``gpu_workers`` (a semaphore over calls would
    not -- contexts would leak to ``threads``). ``make_krill_aligner`` auto-falls
    back to CPU per worker on GPU-init/OOM.
  * Each worker builds ONE runner in the pool initializer and reuses it across every
    interval it is handed; interval indices (not objects) cross the boundary, and
    each worker rebuilds the deterministic interval list itself.
"""
from __future__ import annotations

import atexit
import copy
import logging
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Optional, Tuple

from fin.analysis.abundance_refit import ResponsibilityLedger
from fin.analysis.quantification import QuantResult
from fin.pipeline.config import PipelineConfig

logger = logging.getLogger(__name__)

# Per-worker globals, populated once by _worker_init in each spawned process.
_WORKER_RUNNER = None
_WORKER_INTERVALS: List = []
_WORKER_ID = -1


def _worker_init(config: PipelineConfig, counter, gpu_workers: int, log_level: str) -> None:
    """Pool initializer: claim a worker id, build a runner, rebuild intervals."""
    global _WORKER_RUNNER, _WORKER_INTERVALS, _WORKER_ID

    from fin.io.interval_manager import generate_isolated_intervals
    from fin.pipeline.runner import PipelineRunner
    from fin.utils.log_config import setup_logger

    # Spawned processes start with no logging handlers.
    setup_logger("fin", level=log_level)

    with counter.get_lock():
        _WORKER_ID = counter.value
        counter.value += 1

    role_gpu = _WORKER_ID < gpu_workers
    cfg = copy.deepcopy(config)
    cfg.use_gpu = role_gpu

    _WORKER_RUNNER = PipelineRunner(cfg)
    _WORKER_RUNNER.setup()
    # Spawned workers exit without unwinding run_parallel's scope; release the
    # runner's signal/BAM handles (and any GPU context) at process teardown.
    atexit.register(_WORKER_RUNNER.cleanup)

    result = generate_isolated_intervals(
        cfg.bam_path,
        gtf_path=cfg.gtf_path,
        max_gap=cfg.max_gap,
        max_reads=cfg.max_reads,
    )
    _WORKER_INTERVALS = result["intervals"]
    logger.info("worker %d ready (gpu=%s)", _WORKER_ID, role_gpu)


def _worker_process_index(
    i: int,
) -> Optional[Tuple[List[QuantResult], Optional[ResponsibilityLedger]]]:
    """Process the interval at deterministic index ``i`` in this worker."""
    interval = _WORKER_INTERVALS[i]
    try:
        return _WORKER_RUNNER.process_interval(interval)
    except Exception:
        logger.exception(
            "worker %d failed on interval %d (%s)", _WORKER_ID, i, interval.region_string
        )
        raise


def run_parallel(
    config: PipelineConfig,
    intervals: List,
    threads: int,
    gpu_workers: int,
    log_level: str,
) -> List[Tuple[List[QuantResult], Optional[ResponsibilityLedger]]]:
    """Run ``process_interval`` over all intervals in a spawn process pool.

    Returns non-empty interval ``(results, ledger)`` pairs in completion order.
    Aggregation is mathematically commutative but NOT bitwise associative, so
    completion order perturbs floating-point sums. The refit-enabled path must
    therefore canonicalize this list before aggregating; see
    ``runner._order_interval_outputs``. Do not assume this output is ordered.
    """
    n = len(intervals)
    ctx = mp.get_context("spawn")
    counter = ctx.Value("i", 0)
    all_results: List[
        Tuple[List[QuantResult], Optional[ResponsibilityLedger]]
    ] = []

    ex = ProcessPoolExecutor(
        max_workers=threads,
        mp_context=ctx,
        initializer=_worker_init,
        initargs=(config, counter, gpu_workers, log_level),
    )
    try:
        futures = {ex.submit(_worker_process_index, i): i for i in range(n)}
        done = 0
        for fut in as_completed(futures):
            quant = fut.result()  # re-raises worker exceptions in the parent
            done += 1
            i = futures[fut]
            logger.info(
                "Completed interval %d/%d: %s", done, n, intervals[i].region_string
            )
            if quant is not None:
                all_results.append(quant)
    except BaseException:
        # Fail fast on any worker crash or KeyboardInterrupt: cancel not-yet-started
        # futures and return without blocking on the remaining workload (a plain
        # finally: shutdown(wait=True) would stall the failure behind every queued
        # interval). Already-running tasks cannot be force-killed via the executor;
        # they finish in their own process but the parent no longer waits on them.
        ex.shutdown(wait=False, cancel_futures=True)
        raise
    else:
        ex.shutdown(wait=True)
    return all_results
