# archaea1338__OG0001360 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 121
- Single: 85
- Coexist (2+): 36
- Coexist fraction: 0.297521
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 17
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.84065044, logL=-73.56740833
- Free: r=26.85488, A=0.14293841, logL=-55.06732885
- Delta logL: 18.500079
- Boundary LRT p: 5.90598e-10
- Profile-grid interval: 13.154848786908648 to 75.18211353058093
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 155.732170
- S -> C: 9.809773
- C -> S: 11.041639
- C -> C: 63.416418

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
