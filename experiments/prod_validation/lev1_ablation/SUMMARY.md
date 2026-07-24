# Lever-1 containment ablation (SIRV4, honest metrics, CPU)

Per (sample, mode): off vs on. nout / honest Sn@3 / Pr@3 / F1@3 (and F1@1).

| sample | mode | arm | nout | expr3 | Sn@3 | Pr@3 | F1@3 | F1@1 |
|---|---|---|---|---|---|---|---|---|
| H9_replicate2_run1 | denovo | off | 97 | 81 | 92.6 | 77.3 | 84.3 | 84.2 |
| H9_replicate2_run1 | denovo | on | 96 | 81 | 91.4 | 77.1 | 83.6 | 83.5 |
| H9_replicate2_run1 | guided | off | 97 | 81 | 100.0 | 83.5 | 91.0 | 94.0 |
| H9_replicate2_run1 | guided | on | 96 | 81 | 100.0 | 84.4 | 91.5 | 94.5 |
| H9_replicate2_run2 | denovo | off | 73 | 77 | 75.3 | 79.5 | 77.3 | 74.7 |
| H9_replicate2_run2 | denovo | on | 71 | 77 | 74.0 | 80.3 | 77.0 | 74.4 |
| H9_replicate2_run2 | guided | off | 87 | 77 | 98.7 | 87.4 | 92.7 | 96.5 |
| H9_replicate2_run2 | guided | on | 87 | 77 | 98.7 | 87.4 | 92.7 | 96.5 |
| H9_replicate3_run1 | denovo | off | 93 | 82 | 81.7 | 72.0 | 76.6 | 77.8 |
| H9_replicate3_run1 | denovo | on | 92 | 82 | 81.7 | 72.8 | 77.0 | 78.2 |
| H9_replicate3_run1 | guided | off | 98 | 82 | 98.8 | 82.7 | 90.0 | 93.0 |
| H9_replicate3_run1 | guided | on | 98 | 82 | 98.8 | 82.7 | 90.0 | 93.0 |
| H9_replicate3_run2 | denovo | off | 81 | 76 | 82.9 | 77.8 | 80.3 | 81.0 |
| H9_replicate3_run2 | denovo | on | 78 | 76 | 81.6 | 79.5 | 80.5 | 81.2 |
| H9_replicate3_run2 | guided | off | 92 | 76 | 98.7 | 81.5 | 89.3 | 93.1 |
| H9_replicate3_run2 | guided | on | 91 | 76 | 98.7 | 82.4 | 89.8 | 93.6 |
| H9_replicate4_run1 | denovo | off | 96 | 84 | 85.7 | 75.0 | 80.0 | 76.3 |
| H9_replicate4_run1 | denovo | on | 91 | 84 | 84.5 | 78.0 | 81.1 | 77.2 |
| H9_replicate4_run1 | guided | off | 98 | 84 | 97.6 | 83.7 | 90.1 | 92.9 |
| H9_replicate4_run1 | guided | on | 97 | 84 | 97.6 | 84.5 | 90.6 | 93.3 |
| H9_replicate4_run2 | denovo | off | 69 | 73 | 71.2 | 75.4 | 73.2 | 75.2 |
| H9_replicate4_run2 | denovo | on | 67 | 73 | 69.9 | 76.1 | 72.9 | 74.8 |
| H9_replicate4_run2 | guided | off | 85 | 73 | 97.3 | 83.5 | 89.9 | 94.5 |
| H9_replicate4_run2 | guided | on | 85 | 73 | 97.3 | 83.5 | 89.9 | 94.5 |

## Delta (on - off) per (sample, mode)

| sample | mode | dNout | dSn@3 | dPr@3 | dF1@3 |
|---|---|---|---|---|---|
| H9_replicate2_run1 | denovo | -1 | -1.2 | -0.2 | -0.7 |
| H9_replicate2_run1 | guided | -1 | +0.0 | +0.9 | +0.5 |
| H9_replicate2_run2 | denovo | -2 | -1.3 | +0.8 | -0.3 |
| H9_replicate2_run2 | guided | 0 | +0.0 | +0.0 | +0.0 |
| H9_replicate3_run1 | denovo | -1 | +0.0 | +0.8 | +0.4 |
| H9_replicate3_run1 | guided | 0 | +0.0 | +0.0 | +0.0 |
| H9_replicate3_run2 | denovo | -3 | -1.3 | +1.7 | +0.3 |
| H9_replicate3_run2 | guided | -1 | +0.0 | +0.9 | +0.5 |
| H9_replicate4_run1 | denovo | -5 | -1.2 | +3.0 | +1.1 |
| H9_replicate4_run1 | guided | -1 | +0.0 | +0.9 | +0.5 |
| H9_replicate4_run2 | denovo | -2 | -1.4 | +0.8 | -0.4 |
| H9_replicate4_run2 | guided | 0 | +0.0 | +0.0 | +0.0 |

## Mean over samples

| mode | arm | mean nout | mean Sn@3 | mean Pr@3 | mean F1@3 |
|---|---|---|---|---|---|
| denovo | off | 85 | 81.6 | 76.2 | 78.6 |
| denovo | on | 82 | 80.5 | 77.3 | 78.7 |
| guided | off | 93 | 98.5 | 83.7 | 90.5 |
| guided | on | 92 | 98.5 | 84.1 | 90.8 |
