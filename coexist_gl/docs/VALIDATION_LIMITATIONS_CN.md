# coexist-gl 验证、复现与解释边界

## 1. 数值验证

| Gate | 内容 | 状态 |
|---|---|---|
| 转移核 | `r=1`解析矩阵等于`expm(Q)` | PASS |
| Likelihood | toy pruning等于direct enumeration | PASS |
| Posterior | 事件期望和等于tree edge数 | PASS |
| 聚集检验 | 已知聚集toy得到预期结果 | PASS |
| SpeciesCut | state-blind且抽象likelihood有限 | PASS |

## 2. Actin输入冻结与回归

当前示例为461-tip rooted binary species tree、568-tip archaeal gene tree和568条mapping；354个single、107个coexist，且所有coexist实际均为2 copies。16个gene-tree真核context tips已在冻结前显式排除。

软件0.1.0回归基准：

```text
n_hosts                  461
n_single                 354
n_coexist                107
p_cluster                0.0001
r_hat                    60.7906049076
A_hat                    0.04314541045
logL_null                -249.0791770801
logL_free                -118.2226119172
delta_logL               130.8565651629
p_persistence_one_sided  3.6313597034e-59
```

这些值是软件回归基准，不自动等于论文定稿值。

## 3. 主要限制

### Presence conditioning

零copy hosts被排除。结论条件于“host已有该family”，不估计family在完整物种树上的presence/absence gain/loss。

### Unit edges

所有边长度为1。`A,r`没有每百万年或每替换单位的解释。

### 二状态压缩

所有`n>=2`编码为C。actin只有1或2 copies，其他family可能因3+ copies丢失信息。

### D/T不可分

`G=D+T`。gene tree不进入reconciliation likelihood，所以不能判断coexist由duplication还是transfer产生。

### 边界LRT与profile

单边p值依赖渐近边界分布；小样本建议未来做parametric bootstrap。r区间是离散profile-grid区间，没有端点插值。

### 固定树与mapping

不传播species-tree、gene-tree、mapping和annotation不确定性。copy数还可能受漏注释、contig fragmentation、bin contamination、strain混合和重复预测影响。

## 4. 与ALE关系

本程序计算S/C状态树likelihood，ALE计算gene-tree/species-tree reconciliation likelihood。二者观测空间和归一化不同，absolute logL不可比较；ALE的D/T/L期望也不能与本程序宏观S/C边期望逐项等同。ALE可作外部复杂历史背景，但不是必需输入。

## 5. Cuts与多family

GeneCut不进入主likelihood。SpeciesCut必须state-blind发现；不同cut改变观测空间，absolute likelihood不得跨cut比较。actin当前只有`full` cut，尚无多cut稳健性证据。

批量families需分别对聚集p值和persistence p值做BH-FDR；观测存在不是p值。

## 6. 后续验证优先级

1. alternative species trees和姐妹MAG敏感性；
2. 输出Fitch置换分布分位数和效应量；
3. profile端点连续求根；
4. 小样本parametric bootstrap；
5. 对3+ copy families比较更丰富状态模型；
6. copy state与MAG completeness/contamination关联审计；
7. 批量family FDR。
