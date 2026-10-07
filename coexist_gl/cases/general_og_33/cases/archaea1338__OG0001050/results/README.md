# archaea1338__OG0001050 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 172
- Single: 127
- Coexist (2+): 45
- Coexist fraction: 0.261628
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 33
- Permutations: 9999
- One-sided clustering p: 0.0026
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.7041727, logL=-98.73921863
- Free: r=7.9808677, A=0.28475227, logL=-93.54774413
- Delta logL: 5.1914745
- Boundary LRT p: 0.000635922
- Profile-grid interval: 3.571656212372286 to 16.257276927979312
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 222.979621
- S -> C: 28.172680
- C -> S: 28.662245
- C -> C: 62.185454

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
