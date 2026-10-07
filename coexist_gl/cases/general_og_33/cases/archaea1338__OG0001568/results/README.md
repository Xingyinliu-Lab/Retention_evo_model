# archaea1338__OG0001568 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 167
- Single: 140
- Coexist (2+): 27
- Coexist fraction: 0.161677
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 18
- Permutations: 9999
- One-sided clustering p: 0.0002
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.37961329, logL=-73.64072172
- Free: r=16.777345, A=0.10816548, logL=-61.96698503
- Delta logL: 11.673737
- Boundary LRT p: 6.76105e-07
- Profile-grid interval: 8.143595775458257 to 46.755178074001925
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 266.164332
- S -> C: 12.612239
- C -> S: 12.228730
- C -> C: 40.994698

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
