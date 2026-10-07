# archaea1338__OG0001565 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 130
- Single: 85
- Coexist (2+): 45
- Coexist fraction: 0.346154
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 30
- Permutations: 9999
- One-sided clustering p: 0.0031
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=1.0572794, logL=-83.83682437
- Free: r=10.13666, A=0.35116722, logL=-80.56564554
- Delta logL: 3.2711788
- Boundary LRT p: 0.00526678
- Profile-grid interval: 3.348717643619565 to 22.575329851463824
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 149.494305
- S -> C: 23.185952
- C -> S: 21.147605
- C -> C: 64.172138

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
