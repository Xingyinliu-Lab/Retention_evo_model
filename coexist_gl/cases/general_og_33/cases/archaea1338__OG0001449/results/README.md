# archaea1338__OG0001449 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 134
- Single: 86
- Coexist (2+): 48
- Coexist fraction: 0.358209
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 27
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=1.1122833, logL=-87.37818329
- Free: r=15.010922, A=0.32117826, logL=-79.25911929
- Delta logL: 8.119064
- Boundary LRT p: 2.79293e-05
- Profile-grid interval: 6.359573367528893 to 56.96904293877614
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 163.963556
- S -> C: 23.488701
- C -> S: 14.801409
- C -> C: 63.746335

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
