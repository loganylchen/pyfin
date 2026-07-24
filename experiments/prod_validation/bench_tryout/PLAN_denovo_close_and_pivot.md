# PLAN — close the de-novo precision question, pivot signal-validation to corruption (pyfin's niche)

## Verdict (two independent analyses agree)
De-novo p00 precision parity with isoquant (68→91) is NOT achievable under no-snap.
Oracle math: perfect `c`-truncation removal + perfect fake-junction gate caps at ~84%,
still short of 91%. The residual is real-junction recombination + terminal-extent
ambiguity that neither signal-NLL nor junction adjacency can see. Root blocker = 5'
terminal identifiability in 3'-anchored dRNA (73% of p00 truncations have NO full-length
counterpart in the data — never sequenced). pyfin's defensible win is corruption
robustness (c_jitter honF1 ~31 vs isoquant ~13), from refusing annotation snap.

## Measured/predicted lever ceilings (honest-F1 p00, OFF=28.4)
- graph assembly: -3.0 (over-merge) — TESTED, default-off
- 5'-TSS brake: +0.3 — TESTED, default-off
- containment collapse (Lever-1): +0.5 Pr — TESTED, default-off
- 5'-guarded truncation drop (A): ceiling +1.2 Pr (only 90/391 'c' droppable); Codex +0.5
- signal fake-junction gate (B) on p00: <+1 (only 151/5102); it is a CORRUPTION lever

## Phase 0 — ONE oracle upper-bound study (cheap, decisive gate) [~30 min]
Confirm the ceiling with data before closing. Post-process p00 pyfin GTF (no engine):
1. Oracle-prune: drop every non-'=' multi-exon candidate that is a 3'-suffix sub-chain
   of a longer candidate (the MOST a terminal/containment lever could ever remove),
   keep all '='. Re-run gffcompare + eval_honest on the pruned GTF.
2. Report oracle honest-F1 vs OFF 28.4, and the realistic (5'-TSS-guarded, frozen
   post-EM) subset delta.
GATE: if oracle honF1 gain < +3 → formally CLOSE p00 parity as a goal (record numbers).
If ≥ +3 → the censored-terminal model is worth a real build (Codex predicts +1..+3, so
this is unlikely; proceed only on evidence).

## Phase 1 — pivot signal-validation (lever B) to CORRUPTION [the real payoff]
Codex's key insight: the fake-junction signal gate strengthens pyfin's WINNING niche.
1. Implement `_apply_novel_signal_validation_gate` (pre-EM, novel+multi-exon only) reusing
   the eventalign path (recon map: `read_cand_mean_nll`, `class_junction_window_set`,
   `make_krill_aligner` in m2_junction_nll.py; insertion after junction-dominance gate in
   runner.py process_interval). Default-OFF flag `--novel-signal-gate`.
   Logic: a novel junction whose reads' eventalign NLL fails vs the genome null is a
   misalignment artifact → drop; a well-supported novel junction (real novel splice) → keep.
   This distinguishes real-novel from misaligned — the ONLY tool that can.
2. Measure on CORRUPTION conditions (c_jitter, c_spurious, c_merge) + full + p00/p10.
   Watch the read-reassignment trap (drop few/targeted; prefer frozen post-EM pruning if
   pre-EM regresses).
GATE: honF1 UP on corruption (esp. c_jitter, must stay > isoquant), neutral on full/de-novo,
zero recall loss where possible. Codex-review before any commit.

## Phase 2 — consolidate + reposition
1. Decide commit disposition of the tested default-off experimental code (junction_graph.py
   + 5'-TSS brake + config/cli/runner wiring, 35 tests): commit as documented default-off
   experiment, OR leave uncommitted. All default-off, byte-identical when off.
2. Update pyfin positioning: corruption-robust reference-based dRNA assembly+quant; de-novo
   mode useful but not isoquant-precision-competitive (documented, with the ceiling math).
3. Every commit through the Codex review gate; defaults stay OFF until a corruption-condition
   campaign shows a stable win.

## Non-goals (explicitly dropped, with evidence)
- Chasing p00 struct-precision toward 91% via generation-layer merging/dropping.
- Any annotation/consensus snap (kills the corruption win).
