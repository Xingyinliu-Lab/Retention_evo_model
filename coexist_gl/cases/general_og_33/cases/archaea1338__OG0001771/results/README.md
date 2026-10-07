# archaea1338__OG0001771 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 142
- Single: 102
- Coexist (2+): 40
- Coexist fraction: 0.28169
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 18
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.77905843, logL=-84.32767122
- Free: r=36.581546, A=0.15614768, logL=-63.53989087
- Delta logL: 20.78778
- Boundary LRT p: 5.6701e-11
- Profile-grid interval: 16.24537132784999 to 128.86316854888852
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 200.924150
- S -> C: 13.802413
- C -> S: 7.458632
- C -> C: 59.814805

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
