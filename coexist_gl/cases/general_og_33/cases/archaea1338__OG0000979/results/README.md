# archaea1338__OG0000979 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 146
- Single: 73
- Coexist (2+): 73
- Coexist fraction: 0.5
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 36
- Permutations: 9999
- One-sided clustering p: 0.0003
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=1.9986801, logL=-101.1933629
- Free: r=13.45859, A=0.48978719, logL=-94.24798235
- Delta logL: 6.9453805
- Boundary LRT p: 9.68743e-05
- Profile-grid interval: 6.129539516649347 to 26.855653595522462
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 123.873947
- S -> C: 27.248904
- C -> S: 23.899118
- C -> C: 114.978031

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
