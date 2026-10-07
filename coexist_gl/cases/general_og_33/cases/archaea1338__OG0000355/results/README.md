# archaea1338__OG0000355 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 248
- Single: 194
- Coexist (2+): 54
- Coexist fraction: 0.217742
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 17
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.54749972, logL=-129.5603095
- Free: r=77.682354, A=0.072619864, logL=-72.32043327
- Delta logL: 57.239876
- Boundary LRT p: 5.11528e-27
- Profile-grid interval: 35.25503752861261 to 228.82965272591002
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 385.764648
- S -> C: 12.224248
- C -> S: 6.229987
- C -> C: 89.781116

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
