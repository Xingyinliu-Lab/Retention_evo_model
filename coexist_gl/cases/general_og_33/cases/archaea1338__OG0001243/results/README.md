# archaea1338__OG0001243 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 164
- Single: 114
- Coexist (2+): 50
- Coexist fraction: 0.304878
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 35
- Permutations: 9999
- One-sided clustering p: 0.0023
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.87342326, logL=-100.7815176
- Free: r=8.5250303, A=0.32676631, logL=-96.58693677
- Delta logL: 4.1945808
- Boundary LRT p: 0.00188732
- Profile-grid interval: 3.417464236702818 to 17.95578304573202
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 196.474765
- S -> C: 28.630921
- C -> S: 29.112608
- C -> C: 71.781706

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
