# Lever-1 containment ablation (heya8 dense-locus, honest metrics, CPU)

| sample | mode | arm | nout | expr3 | Sn@3 | Pr@3 | F1@3 | F1@1 |
|---|---|---|---|---|---|---|---|---|
| HEYA8_replicate1_run1 | denovo | off | 80 | 79 | 70.9 | 70.0 | 70.4 | 63.0 |
| HEYA8_replicate1_run1 | denovo | on | 77 | 79 | 69.6 | 71.4 | 70.5 | 63.0 |
| HEYA8_replicate1_run1 | guided | off | 100 | 79 | 98.7 | 78.0 | 87.2 | 94.1 |
| HEYA8_replicate1_run1 | guided | on | 102 | 79 | 98.7 | 76.5 | 86.2 | 94.2 |
| HEYA8_replicate1_run2 | denovo | off | 52 | 53 | 81.1 | 82.7 | 81.9 | 81.0 |
| HEYA8_replicate1_run2 | denovo | on | 50 | 53 | 79.2 | 84.0 | 81.6 | 80.7 |
| HEYA8_replicate1_run2 | guided | off | 73 | 53 | 96.2 | 69.9 | 81.0 | 89.1 |
| HEYA8_replicate1_run2 | guided | on | 73 | 53 | 96.2 | 69.9 | 81.0 | 89.1 |
| HEYA8_replicate2_run1 | denovo | off | 80 | 80 | 65.0 | 65.0 | 65.0 | 56.2 |
| HEYA8_replicate2_run1 | denovo | on | 78 | 80 | 63.8 | 65.4 | 64.6 | 55.7 |
| HEYA8_replicate2_run1 | guided | off | 100 | 80 | 100.0 | 80.0 | 88.9 | 93.7 |
| HEYA8_replicate2_run1 | guided | on | 100 | 80 | 100.0 | 80.0 | 88.9 | 93.7 |
| HEYA8_replicate2_run2 | denovo | off | 63 | 66 | 68.2 | 71.4 | 69.8 | 67.6 |
| HEYA8_replicate2_run2 | denovo | on | 60 | 66 | 66.7 | 73.3 | 69.8 | 67.6 |
| HEYA8_replicate2_run2 | guided | off | 89 | 66 | 97.0 | 71.9 | 82.6 | 91.2 |
| HEYA8_replicate2_run2 | guided | on | 89 | 66 | 97.0 | 71.9 | 82.6 | 91.2 |
| HEYA8_replicate3_run1 | denovo | off | 76 | 78 | 71.8 | 73.7 | 72.7 | 68.7 |
| HEYA8_replicate3_run1 | denovo | on | 72 | 78 | 69.2 | 75.0 | 72.0 | 67.9 |
| HEYA8_replicate3_run1 | guided | off | 98 | 78 | 98.7 | 78.6 | 87.5 | 94.7 |
| HEYA8_replicate3_run1 | guided | on | 99 | 78 | 98.7 | 77.8 | 87.0 | 94.2 |

## Delta (on - off)

| sample | mode | dNout | dSn@3 | dPr@3 | dF1@3 |
|---|---|---|---|---|---|
| HEYA8_replicate1_run1 | denovo | -3 | -1.3 | +1.4 | +0.1 |
| HEYA8_replicate1_run1 | guided | 2 | +0.0 | -1.5 | -1.0 |
| HEYA8_replicate1_run2 | denovo | -2 | -1.9 | +1.3 | -0.4 |
| HEYA8_replicate1_run2 | guided | 0 | +0.0 | +0.0 | +0.0 |
| HEYA8_replicate2_run1 | denovo | -2 | -1.2 | +0.4 | -0.4 |
| HEYA8_replicate2_run1 | guided | 0 | +0.0 | +0.0 | +0.0 |
| HEYA8_replicate2_run2 | denovo | -3 | -1.5 | +1.9 | +0.1 |
| HEYA8_replicate2_run2 | guided | 0 | +0.0 | +0.0 | +0.0 |
| HEYA8_replicate3_run1 | denovo | -4 | -2.6 | +1.3 | -0.7 |
| HEYA8_replicate3_run1 | guided | 1 | +0.0 | -0.8 | -0.5 |

## Mean over samples

| mode | arm | mean nout | mean Sn@3 | mean Pr@3 | mean F1@3 |
|---|---|---|---|---|---|
| denovo | off | 70 | 71.4 | 72.6 | 72.0 |
| denovo | on | 67 | 69.7 | 73.8 | 71.7 |
| guided | off | 92 | 98.1 | 75.7 | 85.4 |
| guided | on | 93 | 98.1 | 75.2 | 85.1 |
