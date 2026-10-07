# archaea1338__OG0001589 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 158
- Single: 124
- Coexist (2+): 34
- Coexist fraction: 0.21519
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 28
- Permutations: 9999
- One-sided clustering p: 0.0669
- Status: `PHYLOGENETIC_CLUSTERING_NOT_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.54847579, logL=-82.27906259
- Free: r=5.7002059, A=0.2759245, logL=-81.3968075
- Delta logL: 0.88225509
- Boundary LRT p: 0.0920315
- Profile-grid interval: 1.0 to 14.174448875753807
- Status: `COEXIST_OBSERVED_NULL_BOUNDARY`

## Posterior宏观边转移期望

- S -> S: 220.794985
- S -> C: 26.954188
- C -> S: 26.197992
- C -> C: 40.052835

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
