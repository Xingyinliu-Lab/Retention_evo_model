# archaea1338__OG0000843 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 246
- Single: 221
- Coexist (2+): 25
- Coexist fraction: 0.101626
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 15
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.21850016, logL=-80.24705563
- Free: r=23.856932, A=0.057891786, logL=-60.93303194
- Delta logL: 19.314024
- Boundary LRT p: 2.56384e-10
- Profile-grid interval: 10.78186348109763 to 80.17320519922865
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 433.867704
- S -> C: 10.898487
- C -> S: 8.559032
- C -> C: 36.674777

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
