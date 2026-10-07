# archaea1338__OG0001103 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 126
- Single: 74
- Coexist (2+): 52
- Coexist fraction: 0.412698
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 27
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=1.4011675, logL=-85.37534755
- Free: r=16.283224, A=0.35033601, logL=-75.75419081
- Delta logL: 9.6211567
- Boundary LRT p: 5.75665e-06
- Profile-grid interval: 7.511892867961312 to 42.72732975153435
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 135.721964
- S -> C: 21.144448
- C -> S: 15.869056
- C -> C: 77.264532

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
