# archaea1338__OG0002388 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 127
- Single: 101
- Coexist (2+): 26
- Coexist fraction: 0.204724
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 22
- Permutations: 9999
- One-sided clustering p: 0.1546
- Status: `PHYLOGENETIC_CLUSTERING_NOT_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.512882, logL=-64.3319039
- Free: r=9.7581562, A=0.18586089, logL=-62.07673157
- Delta logL: 2.2551723
- Boundary LRT p: 0.0168452
- Profile-grid interval: 2.0026477423610856 to 20.705830231467868
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 178.656486
- S -> C: 14.651751
- C -> S: 18.019540
- C -> C: 40.672223

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
