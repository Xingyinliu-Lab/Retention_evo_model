# archaea1338__OG0002523 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 121
- Single: 92
- Coexist (2+): 29
- Coexist fraction: 0.239669
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 18
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.6272599, logL=-66.57542076
- Free: r=18.115498, A=0.17748082, logL=-58.39148112
- Delta logL: 8.1839396
- Boundary LRT p: 2.60809e-05
- Profile-grid interval: 7.27377049069667 to 56.7931968649636
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 176.442183
- S -> C: 13.825040
- C -> S: 9.697516
- C -> C: 40.035261

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
