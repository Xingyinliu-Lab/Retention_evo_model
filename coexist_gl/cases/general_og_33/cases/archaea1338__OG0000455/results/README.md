# archaea1338__OG0000455 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 295
- Single: 237
- Coexist (2+): 58
- Coexist fraction: 0.19661
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 35
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.48321485, logL=-145.8674518
- Free: r=18.137013, A=0.12551819, logL=-120.8761331
- Delta logL: 24.991319
- Boundary LRT p: 7.75562e-13
- Profile-grid interval: 11.133749286154089 to 37.150285547472926
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 448.152834
- S -> C: 24.641576
- C -> S: 24.248364
- C -> C: 90.957225

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
