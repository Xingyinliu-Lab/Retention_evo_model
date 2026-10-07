# archaea1338__OG0000345 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 338
- Single: 299
- Coexist (2+): 39
- Coexist fraction: 0.115385
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 32
- Permutations: 9999
- One-sided clustering p: 0.0013
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.25766295, logL=-120.5235329
- Free: r=7.6296117, A=0.12143368, logL=-114.3719285
- Delta logL: 6.1516045
- Boundary LRT p: 0.00022609
- Profile-grid interval: 3.7481862723421417 to 14.048900332174082
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 560.846263
- S -> C: 29.778864
- C -> S: 32.500375
- C -> C: 50.874499

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
