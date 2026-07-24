#!/bin/bash
cd /SSD/logan/dev/pyfin
OUT=/SSD/logan/dev/pyfin/.rereview
run_one() {
  c=$1
  case $c in 84241b1|d277b4e) kind="ADDITIVE FEATURE (not wired into production)";; *) kind="BEHAVIOR-PRESERVING REFACTOR (byte-identical bar)";; esac
  subj=$(git log -1 --format=%s "$c")
  {
    echo "Re-review a git commit in the pyfin repo with fresh, skeptical eyes. A PRIOR review by a weaker model (gpt-5.5) marked it CLEAN. Independently re-verify with the stronger model and catch anything it may have missed. Read the embedded diff (git show); use a PTY if the non-PTY shell returns blank stdout, but the diff below is self-contained."
    echo
    echo "Commit: $c — $subj"
    echo "Nature: $kind"
    echo
    case $c in
      84241b1|d277b4e)
        echo "Verify: correctness, DETERMINISM, no shared-state mutation, and that it CANNOT change existing production behaviour (cluster_read_chains + discovery path untouched; new symbols only). For the GTF-attach commit specifically: non-bridging (a GTF never merges two read families), zero-read pool integrity (read_pool/variant_reads never get GTF), GTF-only families for unmatched GTF." ;;
      *)
        echo "Acceptance bar: BYTE-IDENTICAL behaviour. Verify the move/extraction changed NO behaviour: no reordered statements that matter, no changed default/threshold, no altered condition or control flow, preserved test/mock seams (patch.object / module-level patch targets), no import cycle, no lost local variable or import, sequential-survivor / lazy-build / fold semantics preserved. Flag ANY line that could shift output." ;;
    esac
    echo
    echo "Give a clear verdict for THIS commit: CLEAN or ISSUES (list them with line cites). Be concrete and skeptical. Review only; do not modify code."
    echo
    echo "===== git show $c ====="
    git show "$c"
  } > "$OUT/prompt_$c.txt"
  codex exec -m gpt-5.6-sol -c model_reasoning_effort="xhigh" --dangerously-bypass-approvals-and-sandbox \
    < "$OUT/prompt_$c.txt" > "$OUT/out_$c.txt" 2>&1
  echo "done $c ($(date +%H:%M:%S))" >> "$OUT/progress.log"
}
export -f run_one
export OUT
printf '%s\n' bc39ce8 b01bcf6 6d99c5b e1698ae cda1a5b 1f75400 84241b1 d277b4e \
  | xargs -P4 -I{} bash -c 'run_one "$@"' _ {}
echo "ALL DONE ($(date +%H:%M:%S))" >> "$OUT/progress.log"
