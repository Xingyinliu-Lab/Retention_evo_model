# archaea1338__OG0000648 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 251
- Single: 210
- Coexist (2+): 41
- Coexist fraction: 0.163347
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 21
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.38237301, logL=-111.2866064
- Free: r=31.552633, A=0.072706494, logL=-78.88399848
- Delta logL: 32.402608
- Boundary LRT p: 4.134e-16
- Profile-grid interval: 16.45962801493781 to 66.77747970168217
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 407.814930
- S -> C: 12.934820
- C -> S: 11.561724
- C -> C: 67.688525

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
