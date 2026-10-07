# archaea1338__OG0000566 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 272
- Single: 237
- Coexist (2+): 35
- Coexist fraction: 0.128676
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 20
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=0.2874258, logL=-103.8229806
- Free: r=25.09315, A=0.065707034, logL=-77.15202032
- Delta logL: 26.67096
- Boundary LRT p: 1.40131e-13
- Profile-grid interval: 12.780914894715938 to 61.66328716053641
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 460.882597
- S -> C: 13.198758
- C -> S: 12.159732
- C -> C: 55.758913

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
