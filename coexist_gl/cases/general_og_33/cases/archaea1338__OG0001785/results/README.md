# archaea1338__OG0001785 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 120
- Single: 80
- Coexist (2+): 40
- Coexist fraction: 0.333333
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 27
- Permutations: 9999
- One-sided clustering p: 0.0051
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.99701587, logL=-76.35356814
- Free: r=9.5350597, A=0.34355395, logL=-72.74913925
- Delta logL: 3.6044289
- Boundary LRT p: 0.00362723
- Profile-grid interval: 3.300409342111571 to 21.744983534927822
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 138.302661
- S -> C: 21.240972
- C -> S: 20.515651
- C -> C: 57.940716

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
