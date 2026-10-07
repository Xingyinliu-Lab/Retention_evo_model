# archaea1338__OG0000244 coexist-gl结果

本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。

## 观测状态

- Hosts: 318
- Single: 193
- Coexist (2+): 125
- Coexist fraction: 0.393082
- Applicability: `FIT_ELIGIBLE`

## 系统发育聚集

- Fitch steps: 52
- Permutations: 9999
- One-sided clustering p: 0.0001
- Status: `PHYLOGENETIC_CLUSTERING_SUPPORTED`

## Coexistence persistence拟合

- Null: r=1, A=1.2906007, logL=-212.9960037
- Free: r=27.976509, A=0.22836673, logL=-168.5581747
- Delta logL: 44.437829
- Boundary LRT p: 2.10213e-21
- Profile-grid interval: 18.16284665559271 to 48.973073570232174
- Status: `COEXIST_PERSISTENCE_SUPPORTED`

## Posterior宏观边转移期望

- S -> S: 361.450592
- S -> C: 35.997271
- C -> S: 29.729068
- C -> C: 206.823069

## 解释边界

这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。
