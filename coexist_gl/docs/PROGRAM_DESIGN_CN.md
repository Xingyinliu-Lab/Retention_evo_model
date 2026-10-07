# coexist-gl 程序设计说明

## 1. 目标与主流程

程序接收analysis-ready species tree、gene tree和gene-to-host mapping，完成输入闭合、copy-state生成、适用性判断、聚集检验、完整species tree拟合、profile/surface、branch posterior、cut审计和provenance写出。它不调用ALE、不推断taxonomy，也不区分D/T来源。

```text
CLI
  -> load_config
  -> validate_inputs
  -> derive S/C states and applicability
  -> gene-tree split audit
  -> Fitch permutation clustering test
  -> full-tree null/free likelihood fit
  -> r profile and r x A surfaces
  -> branch posterior and macro-edge expectations
  -> optional state-blind SpeciesCut fits
  -> summary and completion hashes
```

## 2. 实际代码结构

```text
src/coexist_gl/
├── __init__.py          # 软件/合同版本
├── cli.py               # validate/run入口
├── config.py            # YAML解析与约束
├── inputs.py            # tree/mapping校验和S/C生成
├── clustering.py        # Fitch与tip permutation
├── model.py             # 转移矩阵、pruning、inside/outside
├── fit.py               # null/free、profile、surface轴
├── gene_cuts.py         # gene-tree内部splits
├── species_cuts.py      # state-blind topology cut
├── abstract_model.py    # SpeciesCut beta-binomial likelihood
├── abstract_fit.py      # SpeciesCut拟合
├── runner.py            # 全流程编排与输出
└── utils.py             # SHA、JSON、TSV
```

## 3. 关键对象

- `FamilyConfig`：输入绝对路径、seed、拟合起点、surface和cut设置。
- `ValidatedInputs`：解析后的trees、mapping、copy counts、S/C states和输入SHA。
- `StateTreeLikelihood`：预编译的postorder/preorder索引、children和tip-state数组。
- `Point`：一个`r,A,log_likelihood`拟合点。
- `SpeciesCut/AbstractNode`：可选抽象树，不进入full-tree主结果。

## 4. 模块责任

### 配置与输入

`config.py`强制v2 schema、rooted binary species tree、unit edges、`single=1/coexist_min=2`和`r>=1`。`inputs.py`要求gene-tree tips等于included genes，mapped hosts等于species-tree tips；不做模糊名称匹配或静默删除。

### 概率模型

`model.py`使用二状态解析矩阵指数。`StateTreeLikelihood.evaluate()`默认只做inside；请求事件期望时才做outside pass。所有计算在log space中完成。

### 拟合

`fit.py`中的null是一维bounded `log A`优化；free以多起点L-BFGS-B优化`log A,log r`；profile固定每个r重新优化A；surface复用同一model逐点计算。

### 聚集检验

`clustering.py`预先构造postorder child索引，每次置换仅替换tip state vector，避免重复解析tree。

### Cuts

`gene_cuts.py`提取主gene tree的全部非平凡内部split，仅作结构审计。`species_cuts.py`只读取拓扑和目标block数，不允许根据S/C状态或likelihood寻找有利cut。

### Runner与输出

`runner.py`负责编排和写文件，不复制概率公式。若适用性门失败，保留输入、状态和GeneCut审计后正常退出，不伪造r或p值。

## 5. 性能设计

固定`A,r`的全树likelihood为`O(E)`。461-tip actin树单次计算为毫秒级；完整流程约十余秒，主要工作包括9,999次Fitch置换、优化、profile和100点surface。当前无需复杂family内并行；批量families可在family层独立并行，每个family使用独立输出目录。

## 6. 复现与完成marker

关键输入和结果记录SHA256。`SUMMARY_COMPLETE.json`只在完整流程末尾生成，绑定contract version、config、species tree、gene tree、genes.tsv和summary SHA。

## 7. 测试合同

1. `r=1`解析转移矩阵等于`expm(Q)`；
2. pruning等于toy direct enumeration；
3. posterior事件期望和等于edge count；
4. 聚集toy得到预期Fitch和p值；
5. SpeciesCut state-blind且抽象likelihood有限。

## 8. 尚未实现

- 0/1/2+模型；
- branch-length/dated模型；
- D/T独立估计和reconciliation histories；
- parametric-bootstrap LRT；
- tree/mapping不确定性传播；
- profile端点插值；
- 批量family FDR汇总器；
- taxonomy下载、名称匹配或reference-tree裁剪。

这些缺口不得通过解释性文字伪装成当前程序能力。
