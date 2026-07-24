# Lever-1 ablation on HUMAN GENCODE (raw gffcompare vs full annotation)

OFF = e268c9b SIF output (== live OFF, byte-identical). ON = live --containment-collapse.
Sn/Pr = transcript-level (diluted by full-annotation denominator; DELTA is the signal). c%/=% = query class-code share.

| sample | mode | arm | nout | Tx_Sn | Tx_Pr | c% | =% |
|---|---|---|---|---|---|---|---|
| H9_replicate2_run2 | full | off | 16683 | 6.2 | 94.1 | 0.3 | 95.6 |
| H9_replicate2_run2 | full | on | 16849 | 6.3 | 94.2 | 0.2 | 95.6 |
| H9_replicate2_run2 | nogtf | off | 5574 | 1.4 | 61.7 | 14.7 | 64.2 |
| H9_replicate2_run2 | nogtf | on | 5498 | 1.4 | 62.3 | 13.4 | 64.8 |
| H9_replicate4_run2 | full | off | 14141 | 5.4 | 95.3 | 0.2 | 96.6 |
| H9_replicate4_run2 | full | on | 14326 | 5.4 | 95.3 | 0.1 | 96.6 |
| H9_replicate4_run2 | nogtf | off | 4584 | 1.2 | 63.9 | 15.4 | 66.3 |
| H9_replicate4_run2 | nogtf | on | 4572 | 1.2 | 64.3 | 14.1 | 66.6 |

## Delta (on - off)

| sample | mode | dNout | dTx_Pr | dc% | d=% |
|---|---|---|---|---|---|
| H9_replicate2_run2 | full | 166 | +0.1 | -0.1 | +0.0 |
| H9_replicate2_run2 | nogtf | -76 | +0.6 | -1.3 | +0.6 |
| H9_replicate4_run2 | full | 185 | +0.0 | -0.0 | -0.0 |
| H9_replicate4_run2 | nogtf | -12 | +0.4 | -1.2 | +0.3 |
