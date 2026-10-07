# archaea1338__OG0000477 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 257
- Single: 198
- Coexist (2+): 59
- Coexist fraction: 0.229572
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 37
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.59041654, logL=-138.2320322
- Free: r=17.028984, A=0.16431735, logL=-119.9861442
- Delta logL: 18.245888
- Boundary LRT p: 7.66559e-10
- Profile-grid interval: 9.421574433824947 to 33.271160017846896
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 375.492384
- S -> C: 27.066625
- C -> S: 22.849493
- C -> C: 86.591498

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
