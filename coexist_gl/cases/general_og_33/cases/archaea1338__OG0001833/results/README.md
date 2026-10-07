# archaea1338__OG0001833 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 105
- Single: 63
- Coexist (2+): 42
- Coexist fraction: 0.4
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 23
- Permutations: 9999
- One-sided clustering p: 0.0004
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=1.3297578, logL=-70.64387526
- Free: r=14.32302, A=0.34464777, logL=-63.78781885
- Delta logL: 6.8560564
- Boundary LRT p: 0.000106538
- Profile-grid interval: 6.259680502436774 to 36.610225352896315
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 110.615849
- S -> C: 16.917710
- C -> S: 15.330179
- C -> C: 65.136261

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
