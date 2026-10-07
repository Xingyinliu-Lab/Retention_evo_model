# archaea1338__OG0000688 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 227
- Single: 174
- Coexist (2+): 53
- Coexist fraction: 0.23348
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 23
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.60225057, logL=-123.1045638
- Free: r=33.616536, A=0.078935013, logL=-82.47279672
- Delta logL: 40.631767
- Boundary LRT p: 9.87691e-20
- Profile-grid interval: 19.81844833211058 to 64.42743269600852
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 325.922918
- S -> C: 11.157276
- C -> S: 15.696101
- C -> C: 99.223705

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
