# archaea1338__OG0002358 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 107
- Single: 63
- Coexist (2+): 44
- Coexist fraction: 0.411215
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 31
- Permutations: 9999
- One-sided clustering p: 0.1468
- Status: `PHYLOGENETIC_CLUSTERING_NOT_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=1.396607, logL=-72.47069256
- Free: r=2.4888721, A=0.91273242, logL=-72.44317332
- Delta logL: 0.027519243
- Boundary LRT p: 0.407259
- Profile-grid interval: 1.0 to 10.374716437208074
- Status: `COEXIST_OBSERVED_NULL_BOUNDARY`

## Posterior宏观边转移期望

- S -> S: 88.782422
- S -> C: 37.380786
- C -> S: 36.618703
- C -> C: 49.218090

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
