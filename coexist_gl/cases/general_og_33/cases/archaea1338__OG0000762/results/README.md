# archaea1338__OG0000762 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 236
- Single: 195
- Coexist (2+): 41
- Coexist fraction: 0.173729
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 27
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.41521187, logL=-108.6773411
- Free: r=18.477862, A=0.093475291, logL=-88.13034232
- Delta logL: 20.546999
- Boundary LRT p: 7.25406e-11
- Profile-grid interval: 10.451853734385365 to 35.01324061463501
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 362.243593
- S -> C: 14.795250
- C -> S: 20.284967
- C -> C: 72.676189

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
