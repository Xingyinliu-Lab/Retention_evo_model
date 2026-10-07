# archaea1338__OG0001580 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 107
- Single: 52
- Coexist (2+): 55
- Coexist fraction: 0.514019
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 25
- Permutations: 9999
- One-sided clustering p: 0.0003
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=2.1142741, logL=-74.12149168
- Free: r=14.00514, A=0.49803127, logL=-69.09159246
- Delta logL: 5.0298992
- Boundary LRT p: 0.000757699
- Profile-grid interval: 5.426281518057307 to 31.50538097845153
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 86.877804
- S -> C: 19.666619
- C -> S: 17.415863
- C -> C: 88.039715

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
