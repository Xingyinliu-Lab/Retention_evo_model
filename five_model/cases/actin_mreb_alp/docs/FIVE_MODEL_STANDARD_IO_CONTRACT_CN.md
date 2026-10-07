# Composition-corrected 五模型分析：标准输入/输出合同

## 1. 适用范围

本文档规定三状态、五模型、composition-corrected、unit-topology 分析的标准数据接口。它是分析核心的输入/输出合同，不包含绘图文件，也不包含已经废弃的旧模型、source-branch-length 分析或方向性模型。

当前正式实例为 MreB/ALP aggregate：

- `M-only`：MAG 中存在 M aggregate，但不存在 ALP aggregate；
- `M+ALP coexist`：同一 MAG 中同时存在 M aggregate 与 ALP aggregate；
- `ALP-only`：MAG 中存在 ALP aggregate，但不存在 M aggregate。

M aggregate 包含 `MreB/Mbl-like`、`MreB-like` 与 `C/ALP-N-related`；ALP aggregate 包含 cALP 与 dALP。状态必须在进入本分析前冻结。本模型不读取蛋白序列、不重新分类 family，也不根据 taxonomy 猜测状态。

## 2. 分层接口

分析分成两个严格分离的层次。

### 2.1 输入构建层

输入构建层负责：

1. 读取冻结 MAG 三状态表；
2. 将 MAG 映射到各物种树叶标签；
3. 对每棵物种树保留具有确定状态的叶；
4. 将替代物种树的操作根投影到 GTDB guide root；
5. 检查树叶集合与状态键集合完全相等；
6. 写出唯一的分析就绪 JSON。

当前构建入口：

```text
code/run_mreb_alp_composition_corrected_five_models.py
```

该入口同时包含输入构建和拟合能力，但正式复算时应优先使用第 2.2 节的冻结输入入口，避免上游映射变化悄然改变分析集。

### 2.2 核心拟合层

核心拟合只允许读取：

```text
input/six_tree_states.json
```

正式冻结复算入口：

```text
code/fit_frozen_input_no_biopython.py
```

核心拟合不得再读取原始 taxonomy、MAG metadata、蛋白表或外部物种树，不进行自动 pruning、名称修复或状态重判。

## 3. 标准输入合同

### 3.1 文件级要求

唯一必需输入是 UTF-8 JSON：

```text
input/six_tree_states.json
```

顶层结构：

```json
{
  "schema_version": "mreb_alp_composition_corrected_five_model_v1",
  "branch_scale": "unit_topology",
  "trees": [ ... ]
}
```

约束：

- `schema_version` 必须精确匹配；
- `branch_scale` 必须为 `unit_topology`；
- `trees` 至少含一棵树，当前正式实例为六棵；
- 每个 `name` 在文件内唯一；
- 不允许缺失状态、重复 tip 或未映射 tip。

### 3.2 每棵树对象

每个 `trees[]` 元素必须包含：

| 字段 | 类型 | 必需 | 定义 |
|---|---:|---:|---|
| `name` | string | 是 | 树的唯一分析 ID |
| `newick` | string | 是 | 已裁剪、已确定操作根的 Newick |
| `states` | object | 是 | `tip_id -> state` 映射 |
| `source_tree` | string | 是 | 原始物种树绝对路径，仅作 provenance |
| `root_semantics` | string | 是 | 当前正式值为 `gtdb_projected` |
| `root_projection_mismatch` | integer | 是 | 操作根投影与 guide split 的 mismatch 计数 |
| `n_tips` | integer | 是 | Newick 叶数 |
| `state_counts` | object | 是 | 三状态计数审计 |

`states` 只允许三个精确字符串：

```text
M-only
M+ALP coexist
ALP-only
```

### 3.3 集合闭合要求

对每棵树必须同时满足：

```text
set(Newick tips) == set(states.keys())
n_tips == number of unique Newick tips
n_tips == number of states entries
sum(state_counts.values()) == n_tips
```

任何一项失败均为硬错误，不得进入优化。

### 3.4 树要求

- 树必须具有一个明确的操作根；这个根用于离散性状似然递归，不自动等同于绝对生物学根。
- tip 名必须唯一且与 `states` 键逐字符一致。
- 内部节点允许二叉或多叉；不得存在空子树或仅含一个后代的冗余内部节点。
- 当前分析忽略 Newick branch lengths：每条保留边统一设为一个离散时间单位。
- 六棵树是拓扑敏感性分析，不是六个独立生物学重复。

### 3.5 当前冻结实例

| tree | tips | M-only | coexist | ALP-only |
|---|---:|---:|---:|---:|
| GTDB_r226 | 4,208 | 2,546 | 191 | 1,471 |
| A1338_addcore | 284 | 64 | 31 | 189 |
| A1338_denovo | 284 | 64 | 31 | 189 |
| A1338_ASTRAL64 | 284 | 64 | 31 | 189 |
| Asgard971_addcore | 851 | 14 | 182 | 655 |
| Asgard971_ASTRAL89 | 851 | 14 | 182 | 655 |

冻结输入：

```text
file: input/six_tree_states.json
bytes: 486548
lines: 6858
SHA256: 036016315e01fc758c9bd05687dcfda24833456243bc7ea84e6ed5ca8f196de1
```

## 4. 五个模型的固定定义

状态编号固定为：

```text
0 = M-only
1 = M+ALP coexist
2 = ALP-only
```

所有速率均为非负；优化在 log-rate 空间进行。

### 4.1 ARD unrestricted

允许全部六个有向转换：

```text
0->1, 1->0, 1->2, 2->1, 0->2, 2->0
```

- 自由参数数 `k=6`；
- root prior 为拟合 Q 的 stationary distribution；
- 这是灵活 benchmark，不自动代表独立生物学机制。

### 4.2 Sequential retention

正式代码名：`adjacency_coexist_middle`。

```text
M-only <-> M+ALP coexist <-> ALP-only
```

- 允许 `0<->1` 与 `1<->2`；
- 禁止 `0<->2` 直接转换；
- `k=4`；
- root prior 为 stationary(Q)；
- 模型可逆，因此不区分 M-first 与 ALP-first 方向。

### 4.3 ALP-middle adjacency

正式代码名：`adjacency_ALP_middle`。

```text
M-only <-> ALP-only <-> M+ALP coexist
```

- 允许 `0<->2` 与 `2<->1`；
- `k=4`；
- root prior 为 stationary(Q)。

### 4.4 M-middle adjacency

正式代码名：`adjacency_M_middle`。

```text
ALP-only <-> M-only <-> M+ALP coexist
```

- 允许 `2<->0` 与 `0<->1`；
- `k=4`；
- root prior 为 stationary(Q)。

### 4.5 Ancestral coexistence strict

正式代码名：`ancestral_coexistence_strict`。

```text
root = M+ALP coexist
M+ALP coexist -> M-only
M+ALP coexist -> ALP-only
```

- 只允许 `1->0` 与 `1->2`；
- 两个 single states 为吸收态；
- `k=2`；
- root prior 固定为 coexist，不从 Q 的 stationary distribution 推导。

## 5. 似然与优化合同

### 5.1 转移矩阵

对每个模型建立三状态连续时间 Markov rate matrix `Q`：

```text
Q[i,j] = fitted rate for an allowed i->j transition
Q[i,j] = 0 for a forbidden transition
Q[i,i] = -sum(Q[i,j], j != i)
```

因为使用 unit topology，每条物种树边共享：

```text
P = expm(Q)
```

原始 Newick branch length 不进入似然。

### 5.2 树似然

采用标准 pruning recursion：tip 为确定的 one-hot 三状态向量；内部节点对子树条件似然逐子节点相乘；每个内部节点做数值缩放并累计 log-scale。根部使用该模型规定的 root prior。

输出：

```text
log_likelihood = log P(observed tip states | tree, fitted Q, root policy)
```

### 5.3 优化

- 优化器：SciPy `L-BFGS-B`；
- 参数：允许边的 log-rate；
- bounds：每个 log-rate 为 `[-10, 10]`；
- 五个固定初始点：全 0、全 -1、全 +1、`linspace(-2,2)`、反向 `linspace(2,-2)`；
- 每个初始点最大 400 iterations；
- 取成功且有限解中的最高似然；若无成功标记，则保留所有尝试中目标值最小者，但 `success` 必须如实为 `FALSE`。

数值优化在不同平台或运行间可有轻微随机/浮点波动，不改变输入合同、模型定义或主要比较。

### 5.4 AIC

```text
AIC = 2*k - 2*log_likelihood
delta_AIC(tree, model) = AIC(tree, model) - min_model AIC(tree, model)
```

其中 `k` 只计自由转换速率；stationary root prior 由 Q 派生，不另计参数；strict 模型的固定 root state 也不计参数。

只允许在同一棵树、同一 tip/state 集合、同一 branch scale 下比较模型 AIC。不得跨树比较 absolute log-likelihood 或 absolute AIC。

## 6. 核心输出合同

### 6.1 主拟合表

路径：

```text
results/MreB_ALP_composition_corrected_five_models_six_trees.tsv
```

粒度：每棵树 × 每个模型一行；当前应为 `6*5=30` 行数据，加一行表头。

| 字段 | 类型 | 定义 |
|---|---:|---|
| `tree` | string | 输入树 ID |
| `model` | string | 五个固定模型代码之一 |
| `root_mode` | string | `stationary` 或固定的 `M+ALP coexist` |
| `branch_scale` | string | 必须为 `unit_topology` |
| `k` | integer | 自由速率参数数 |
| `log_likelihood` | float | 最大化后的树状态 log-likelihood |
| `AIC` | float | `2*k-2*log_likelihood` |
| `success` | boolean string | 优化器成功状态 |
| `gradient_max` | float | 最优解处最大绝对梯度分量，数值诊断 |
| `rates` | string | `state->state:value` 的分号分隔速率 |
| `mapped_MAGs` | integer | 本树参与拟合的 tip/MAG 数 |
| `root_projection_mismatch` | integer | 操作根投影审计量 |
| `n_M_only` | integer | M-only tip 数 |
| `n_coexist` | integer | coexist tip 数 |
| `n_ALP_only` | integer | ALP-only tip 数 |
| `delta_AIC` | float | 同树内相对最佳模型的 AIC 差 |

完成条件：

```text
30/30 rows present
all tree-model keys unique
all log_likelihood and AIC finite
all success == TRUE
all branch_scale == unit_topology
all delta_AIC >= 0 within floating tolerance
one delta_AIC == 0 per tree, allowing numerical ties
```

当前文件冻结信息：

```text
bytes: 9201
lines: 31
SHA256: 0747b9ad2999b5e5244e56768cebd766e6f549b46ea5f6162df7980f9df165f9
```

### 6.2 完成 marker

路径：

```text
FIT_COMPLETE.marker
```

至少包含：

```text
status=complete
models=5
trees=6
branch_scale=unit_topology
```

marker 不能替代对主拟合表的行数、有限值和唯一键审计。

## 7. 标准派生输出

派生入口：

```text
code/analyze_pvalues_and_ard_contributions.py
```

它必须只读取冻结 JSON 与主拟合表，不重新优化五模型。

### 7.1 配对模型比较

路径：

```text
results/paired_tests_vs_sequential_retention.tsv
```

粒度：Sequential retention 与其余四模型各一行。

| 字段 | 定义 |
|---|---|
| `benchmark` | 固定为 `adjacency_coexist_middle` |
| `comparison_model` | 被比较模型 |
| `n_trees` | 成对树数，当前为 6 |
| `retention_better_trees` | AIC(retention) 更低的树数 |
| `comparison_better_trees` | AIC(comparison) 更低的树数 |
| `median_AIC_other_minus_retention` | 六树配对 AIC 差中位数；正值有利于 retention |
| `exact_p_two_sided` | 未校正、双侧、精确配对 signed-rank P 值 |

六棵替代树是敏感性重建，不是独立生物学样本；P 值仅概括跨拓扑一致性，不应表述为基于六个独立重复的总体显著性。

### 7.2 ARD 转换组 posterior moments

路径：

```text
results/MreB_ALP_ARD_transition_group_posterior_moments_six_trees.tsv
```

ARD 的六条转换被分为：

```text
coexist_mediated = 0->1, 1->0, 1->2, 2->1
direct_endpoint  = 0->2, 2->0
```

每棵树输出两行，字段为：

| 字段 | 定义 |
|---|---|
| `tree` | 树 ID |
| `transition_group` | 两个转换组之一 |
| `posterior_mean` | fixed-Q stochastic-history 条件期望转换数 |
| `posterior_variance` | 条件方差 |
| `posterior_sd` | 条件标准差 |
| `log_likelihood_recomputed` | 从冻结 Q 重新计算的 log-likelihood |
| `fit_log_likelihood` | 主拟合表的 ARD log-likelihood |
| `log_likelihood_abs_difference` | 两次似然的绝对差，用于复现审计 |

这不是重新拟合模型，也不是 D/T/L reconciliation。它是在已拟合 ARD Q 下，对离散状态转换路径进行 posterior moment 分解。

### 7.3 ARD 贡献汇总

路径：

```text
results/MreB_ALP_ARD_transition_contribution_summary_six_trees.tsv
```

每棵树一行：

| 字段 | 定义 |
|---|---|
| `coexist_mediated_mean` | 经 coexist 的期望转换数 |
| `direct_endpoint_mean` | M-only 与 ALP-only 直接转换的期望数 |
| `coexist_mediated_fraction` | 前者占两组总量的比例 |
| `direct_endpoint_fraction` | 后者占两组总量的比例 |
| `expected_total` | 两组期望数之和 |
| `expected_total_per_branch` | `expected_total / number_of_retained_edges` |

`coexist_mediated_fraction` 高表示 ARD 拟合的大部分转换也经过 coexist；它不证明转换方向，不等同于 retention 模型的模型选择胜出。

## 8. 非输出与禁止解释

本五模型分析不输出：

- duplication、transfer、loss 的逐事件判定；
- gene-tree/species-tree reconciliation；
- 每条 gene lineage 的历史；
- cALP/dALP 或 MreB/ALP 的祖先—衍生方向；
- duplication 时间；
- HGT donor/recipient；
- 基因功能变化；
- 同一 MAG 多拷贝由 duplication 还是 transfer 产生的确认。

因此不得把 `rates` 或 ARD posterior transitions 改写为 D/T/L 数量，也不得将物种树状态 CTMC 的 likelihood 与 ALE reconciliation 的 absolute likelihood 相减或做 LRT。

## 9. 标准运行顺序

在输入 JSON 已冻结的情况下：

```powershell
python code\fit_frozen_input_no_biopython.py
python code\analyze_pvalues_and_ard_contributions.py
```

标准验收顺序：

1. 核验输入 JSON SHA256、schema、六树名称和集合闭合；
2. 运行五模型拟合；
3. 核验主表 30/30、唯一键、有限值、`success=TRUE`；
4. 核验每棵树的 `delta_AIC` 在该树内部计算；
5. 运行配对检验与 ARD 贡献分解；
6. 核验 ARD 重算似然与主拟合表一致；
7. 最后写/核验 completion marker。

## 10. 正式代码边界

正式核心代码：

```text
code/run_mreb_alp_composition_corrected_five_models.py
code/fit_frozen_input_no_biopython.py
code/analyze_pvalues_and_ard_contributions.py
```

其中：

- `run_mreb...py`：负责从权威上游资产构建冻结 JSON，也可完成拟合；
- `fit_frozen...py`：只基于冻结 JSON 完成正式五模型拟合，是推荐复算入口；
- `analyze...py`：生成配对检验和 ARD 贡献，不改变拟合结果。

绘图脚本和 `figures/`、`source_data/` 不属于本输入/输出合同。
