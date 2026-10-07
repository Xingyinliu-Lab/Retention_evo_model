# archaea1338__OG0001791 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 117
- Single: 85
- Coexist (2+): 32
- Coexist fraction: 0.273504
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 17
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.74496463, logL=-68.5248738
- Free: r=19.681022, A=0.20235428, logL=-57.24955885
- Delta logL: 11.275315
- Boundary LRT p: 1.02339e-06
- Profile-grid interval: 8.046106303541352 to None
- Status: `R_NOT_IDENTIFIABLE`

## Posterior宏观边转移期望

- S -> S: 163.894922
- S -> C: 14.675673
- C -> S: 9.398771
- C -> C: 44.030633

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
