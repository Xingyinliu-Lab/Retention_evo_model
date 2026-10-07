# archaea1338__OG0000710 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 196
- Single: 130
- Coexist (2+): 66
- Coexist fraction: 0.336735
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 39
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=1.0128399, logL=-125.1667966
- Free: r=13.502124, A=0.28672745, logL=-114.1514493
- Delta logL: 11.015347
- Boundary LRT p: 1.34162e-06
- Profile-grid interval: 7.506643249434815 to 25.160150967646466
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 234.255149
- S -> C: 29.487860
- C -> S: 26.969068
- C -> C: 99.287923

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
