# archaea1338__OG0000701 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 195
- Single: 135
- Coexist (2+): 60
- Coexist fraction: 0.307692
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 36
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.88628251, logL=-120.3088142
- Free: r=15.880367, A=0.22332971, logL=-108.0160665
- Delta logL: 12.292748
- Boundary LRT p: 3.55424e-07
- Profile-grid interval: 9.162040463719661 to 29.72012041197959
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 242.570255
- S -> C: 23.961287
- C -> S: 24.707687
- C -> C: 96.760771

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
