# archaea1338__OG0000637 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 197
- Single: 112
- Coexist (2+): 85
- Coexist fraction: 0.431472
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 41
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=1.515096, logL=-134.6679723
- Free: r=20.035707, A=0.33529863, logL=-118.7112205
- Delta logL: 15.956752
- Boundary LRT p: 8.05963e-09
- Profile-grid interval: 10.806720574882696 to 39.584664709196616
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 209.228419
- S -> C: 31.328928
- C -> S: 22.057677
- C -> C: 129.384977

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
