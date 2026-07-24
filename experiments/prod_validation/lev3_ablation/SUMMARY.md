# Lever-3 (mono gate) ablation — off vs on

ON = --drop-mono-exon-novel --min-mono-exon-reads 3 --min-mono-exon-length 200.
mono% = single-exon share of output (Lever-3 target, should drop). sirv/heya8: Sn3/Pr3/F1@3 honest. gencode: Tx_Sn/Pr raw.

| ds | sample | mode | arm | nout | mono% | Sn/Sn3 | Pr/Pr3 | F1@3 |
|---|---|---|---|---|---|---|---|---|
| sirv | H9_replicate2_run1 | denovo | off | 97 | 25.8 | 92.6 | 77.3 | 84.3 |
| sirv | H9_replicate2_run1 | denovo | on | 96 | 25.0 | 91.4 | 77.1 | 83.6 |
| sirv | H9_replicate2_run1 | guided | off | 97 | 33.0 | 100.0 | 83.5 | 91.0 |
| sirv | H9_replicate2_run1 | guided | on | 97 | 33.0 | 100.0 | 83.5 | 91.0 |
| sirv | H9_replicate2_run2 | denovo | off | 73 | 26.0 | 75.3 | 79.5 | 77.3 |
| sirv | H9_replicate2_run2 | denovo | on | 73 | 26.0 | 75.3 | 79.5 | 77.3 |
| sirv | H9_replicate2_run2 | guided | off | 87 | 28.7 | 98.7 | 87.4 | 92.7 |
| sirv | H9_replicate2_run2 | guided | on | 87 | 28.7 | 98.7 | 87.4 | 92.7 |
| sirv | H9_replicate3_run1 | denovo | off | 93 | 25.8 | 81.7 | 72.0 | 76.6 |
| sirv | H9_replicate3_run1 | denovo | on | 92 | 25.0 | 80.5 | 71.7 | 75.9 |
| sirv | H9_replicate3_run1 | guided | off | 98 | 29.6 | 98.8 | 82.7 | 90.0 |
| sirv | H9_replicate3_run1 | guided | on | 98 | 29.6 | 98.8 | 82.7 | 90.0 |
| sirv | H9_replicate3_run2 | denovo | off | 81 | 25.9 | 82.9 | 77.8 | 80.3 |
| sirv | H9_replicate3_run2 | denovo | on | 80 | 25.0 | 81.6 | 77.5 | 79.5 |
| sirv | H9_replicate3_run2 | guided | off | 92 | 30.4 | 98.7 | 81.5 | 89.3 |
| sirv | H9_replicate3_run2 | guided | on | 92 | 30.4 | 98.7 | 81.5 | 89.3 |
| sirv | H9_replicate4_run1 | denovo | off | 96 | 26.0 | 85.7 | 75.0 | 80.0 |
| sirv | H9_replicate4_run1 | denovo | on | 95 | 25.3 | 84.5 | 74.7 | 79.3 |
| sirv | H9_replicate4_run1 | guided | off | 98 | 31.6 | 97.6 | 83.7 | 90.1 |
| sirv | H9_replicate4_run1 | guided | on | 98 | 31.6 | 97.6 | 83.7 | 90.1 |
| sirv | H9_replicate4_run2 | denovo | off | 69 | 27.5 | 71.2 | 75.4 | 73.2 |
| sirv | H9_replicate4_run2 | denovo | on | 69 | 27.5 | 71.2 | 75.4 | 73.2 |
| sirv | H9_replicate4_run2 | guided | off | 85 | 27.1 | 97.3 | 83.5 | 89.9 |
| sirv | H9_replicate4_run2 | guided | on | 85 | 27.1 | 97.3 | 83.5 | 89.9 |
| heya8 | HEYA8_replicate1_run1 | denovo | off | 80 | 8.8 | 70.9 | 70.0 | 70.4 |
| heya8 | HEYA8_replicate1_run1 | denovo | on | 79 | 7.6 | 69.6 | 69.6 | 69.6 |
| heya8 | HEYA8_replicate1_run2 | denovo | off | 52 | 11.5 | 81.1 | 82.7 | 81.9 |
| heya8 | HEYA8_replicate1_run2 | denovo | on | 51 | 9.8 | 79.2 | 82.4 | 80.8 |
| heya8 | HEYA8_replicate2_run1 | denovo | off | 80 | 16.2 | 65.0 | 65.0 | 65.0 |
| heya8 | HEYA8_replicate2_run1 | denovo | on | 79 | 15.2 | 63.8 | 64.6 | 64.2 |
| heya8 | HEYA8_replicate2_run2 | denovo | off | 63 | 15.9 | 68.2 | 71.4 | 69.8 |
| heya8 | HEYA8_replicate2_run2 | denovo | on | 62 | 14.5 | 66.7 | 71.0 | 68.8 |
| heya8 | HEYA8_replicate3_run1 | denovo | off | 76 | 13.2 | 71.8 | 73.7 | 72.7 |
| heya8 | HEYA8_replicate3_run1 | denovo | on | 74 | 10.8 | 69.2 | 73.0 | 71.1 |
| gencode | H9_replicate2_run2 | nogtf | off | 5574 | 12.8 | 1.4 | 61.7 | 0.0 |
| gencode | H9_replicate2_run2 | nogtf | on | 5556 | 12.5 | 1.4 | 61.9 | 0.0 |
| gencode | H9_replicate4_run2 | nogtf | off | 4584 | 12.7 | 1.2 | 63.9 | 0.0 |
| gencode | H9_replicate4_run2 | nogtf | on | 4573 | 12.5 | 1.2 | 64.0 | 0.0 |

## Delta (on - off)

| ds | sample | mode | dNout | dmono% | dSn | dPr | dF1@3 |
|---|---|---|---|---|---|---|---|
| sirv | H9_replicate2_run1 | denovo | -1 | -0.8 | -1.2 | -0.2 | -0.7 |
| sirv | H9_replicate2_run1 | guided | 0 | +0.0 | +0.0 | +0.0 | +0.0 |
| sirv | H9_replicate2_run2 | denovo | 0 | +0.0 | +0.0 | +0.0 | +0.0 |
| sirv | H9_replicate2_run2 | guided | 0 | +0.0 | +0.0 | +0.0 | +0.0 |
| sirv | H9_replicate3_run1 | denovo | -1 | -0.8 | -1.2 | -0.3 | -0.7 |
| sirv | H9_replicate3_run1 | guided | 0 | +0.0 | +0.0 | +0.0 | +0.0 |
| sirv | H9_replicate3_run2 | denovo | -1 | -0.9 | -1.3 | -0.3 | -0.8 |
| sirv | H9_replicate3_run2 | guided | 0 | +0.0 | +0.0 | +0.0 | +0.0 |
| sirv | H9_replicate4_run1 | denovo | -1 | -0.8 | -1.2 | -0.3 | -0.7 |
| sirv | H9_replicate4_run1 | guided | 0 | +0.0 | +0.0 | +0.0 | +0.0 |
| sirv | H9_replicate4_run2 | denovo | 0 | +0.0 | +0.0 | +0.0 | +0.0 |
| sirv | H9_replicate4_run2 | guided | 0 | +0.0 | +0.0 | +0.0 | +0.0 |
| heya8 | HEYA8_replicate1_run1 | denovo | -1 | -1.2 | -1.3 | -0.4 | -0.8 |
| heya8 | HEYA8_replicate1_run2 | denovo | -1 | -1.7 | -1.9 | -0.3 | -1.1 |
| heya8 | HEYA8_replicate2_run1 | denovo | -1 | -1.1 | -1.2 | -0.4 | -0.8 |
| heya8 | HEYA8_replicate2_run2 | denovo | -1 | -1.4 | -1.5 | -0.5 | -1.0 |
| heya8 | HEYA8_replicate3_run1 | denovo | -2 | -2.3 | -2.6 | -0.7 | -1.7 |
| gencode | H9_replicate2_run2 | nogtf | -18 | -0.3 | +0.0 | +0.2 | +0.0 |
| gencode | H9_replicate4_run2 | nogtf | -11 | -0.2 | +0.0 | +0.1 | +0.0 |

## Mean delta by (ds, mode)

| ds | mode | dNout | dmono% | dSn | dPr | dF1@3 |
|---|---|---|---|---|---|---|
| gencode | nogtf | -14 | -0.2 | +0.0 | +0.1 | +0.0 |
| heya8 | denovo | -1 | -1.5 | -1.7 | -0.5 | -1.1 |
| sirv | denovo | -1 | -0.5 | -0.8 | -0.2 | -0.5 |
| sirv | guided | 0 | +0.0 | +0.0 | +0.0 | +0.0 |
