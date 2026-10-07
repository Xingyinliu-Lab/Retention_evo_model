# archaea1338__OG0001088 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 122
- Single: 67
- Coexist (2+): 55
- Coexist fraction: 0.45082
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 24
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=1.6373709, logL=-83.94915069
- Free: r=20.494109, A=0.31226817, logL=-70.19516871
- Delta logL: 13.753982
- Boundary LRT p: 7.82245e-08
- Profile-grid interval: 10.146111000006231 to 46.41771180985006
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 116.269571
- S -> C: 15.879804
- C -> S: 16.144221
- C -> C: 93.706404

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
