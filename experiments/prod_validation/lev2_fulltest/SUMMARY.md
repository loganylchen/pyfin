# Lever-2 FULL test (--novel-junction-min-reads 2), on vs off, all datasets/ratios

sirv/heya8: honest Sn3/Pr3/F1@3 (nanocount expressed-truth). gencode: raw Tx_Sn/Pr.

## Mean delta (on - off) by (dataset, ratio)

| ds | ratio | n | dNout | dSn | dPr | dF1@3 |
|---|---|---|---|---|---|---|
| gencode | c_flip20 | 5 | -141 | +0.0 | +0.7 | +0.0 |
| gencode | c_ir20 | 5 | -155 | +0.0 | +0.8 | +0.0 |
| gencode | c_jitter20_10bp | 5 | -280 | -0.0 | +0.6 | +0.0 |
| gencode | c_merge20 | 4 | -122 | +0.0 | +0.6 | +0.0 |
| gencode | c_skip20 | 5 | -162 | +0.0 | +0.8 | +0.0 |
| gencode | c_spurious20 | 5 | -153 | +0.0 | +0.8 | +0.0 |
| gencode | full | 5 | -148 | +0.0 | +0.8 | +0.0 |
| gencode | nogtf | 5 | -491 | -0.1 | +3.8 | +0.0 |
| gencode | p10 | 5 | -411 | -0.0 | +3.1 | +0.0 |
| gencode | p50 | 5 | -249 | -0.0 | +1.5 | +0.0 |
| gencode | p90 | 4 | -135 | +0.0 | +0.9 | +0.0 |
| gencode | p99 | 3 | -118 | +0.0 | +0.7 | +0.0 |
| heya8 | c_flip5 | 5 | -0 | +0.0 | +0.3 | +0.2 |
| heya8 | c_ir10 | 5 | -1 | -1.1 | +0.2 | -0.3 |
| heya8 | c_jitter10bp | 5 | -1 | +0.0 | +0.6 | +0.4 |
| heya8 | c_merge5 | 5 | -0 | +0.0 | +0.3 | +0.2 |
| heya8 | c_skip10 | 5 | -1 | +0.0 | +0.5 | +0.3 |
| heya8 | c_spurious5 | 5 | -0 | +0.0 | +0.3 | +0.2 |
| heya8 | full | 5 | -0 | +0.0 | +0.3 | +0.2 |
| heya8 | nogtf | 5 | -3 | -1.4 | +1.4 | -0.1 |
| heya8 | p00 | 5 | -3 | -1.4 | +1.4 | -0.1 |
| sirv | c_flip20 | 6 | -1 | +0.0 | +0.5 | +0.3 |
| sirv | c_flip5 | 6 | -1 | +0.0 | +0.6 | +0.3 |
| sirv | c_ir10 | 6 | -1 | +0.0 | +0.6 | +0.3 |
| sirv | c_ir20 | 6 | -2 | -1.3 | +0.1 | -0.4 |
| sirv | c_jitter10bp | 6 | -1 | -0.8 | +0.4 | -0.1 |
| sirv | c_jitter20_10bp | 6 | -2 | -2.1 | +0.2 | -0.8 |
| sirv | c_merge20 | 6 | -1 | +0.0 | +0.6 | +0.4 |
| sirv | c_merge5 | 6 | -2 | -0.8 | +0.8 | +0.1 |
| sirv | c_skip10 | 6 | -1 | +0.0 | +0.6 | +0.3 |
| sirv | c_skip20 | 6 | -1 | -1.3 | +0.1 | -0.5 |
| sirv | c_spurious20 | 6 | -1 | +0.0 | +0.6 | +0.4 |
| sirv | c_spurious5 | 6 | -1 | +0.0 | +0.6 | +0.3 |
| sirv | full | 6 | -1 | +0.0 | +0.6 | +0.3 |
| sirv | p00 | 6 | -3 | -2.1 | +1.0 | -0.4 |
| sirv | p10 | 6 | -3 | -2.1 | +1.0 | -0.4 |
| sirv | p50 | 6 | -2 | -0.8 | +0.6 | -0.0 |
| sirv | p90 | 6 | -2 | -0.8 | +0.9 | +0.2 |
| sirv | p99 | 6 | -1 | +0.0 | +0.6 | +0.3 |

## Mean delta by dataset (all ratios)

| ds | cells | dSn | dPr | dF1@3 |
|---|---|---|---|---|
| sirv | 108 | -0.7 | +0.6 | +0.0 |
| heya8 | 45 | -0.4 | +0.6 | +0.1 |
| gencode | 56 | -0.0 | +1.3 | +0.0 |

## pyfin-ON vs competitors (global mean honest F1@3)

### sirv  (pyfin OFF measured 85.0 vs table pyfin_prod 85.0)

| rank | tool | F1@3 |
|---|---|---|
| 1 | pyfin_ON(lever2) | 85.0 ⬅ |
| 2 | pyfin_OFF(prod) | 85.0 ⬅ |
| 3 | espresso | 78.7 |
| 4 | stringtie3 | 76.6 |
| 5 | isoquant | 72.5 |
| 6 | lafite | 68.8 |
| 7 | flair | 59.3 |
| 8 | bambu | 50.4 |
| 9 | talon | 29.7 |
| 10 | isotools | 27.1 |

### heya8  (pyfin OFF measured 80.2 vs table pyfin_prod ?)

| rank | tool | F1@3 |
|---|---|---|
| 1 | pyfin_ON(lever2) | 80.3 ⬅ |
| 2 | pyfin_OFF(prod) | 80.2 ⬅ |

