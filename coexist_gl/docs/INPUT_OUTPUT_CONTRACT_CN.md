# coexist-gl 输入输出合同

## 1. 必需输入

```text
input/
├── family.yaml
├── species_tree.nwk
├── gene_tree.nwk
└── genes.tsv
```

### `family.yaml`

```yaml
schema_version: coexist_gl_family_v2
family_id: example_family
species_tree:
  path: species_tree.nwk
  rooted: true
  edge_policy: unit
gene_tree:
  path: gene_tree.nwk
  rootedness: unrooted
gene_metadata:
  path: genes.tsv
  include_column: include
copy_state:
  single: 1
  coexist_min: 2
clustering_test:
  permutations: 9999
  seed: 1
fit:
  r_lower_bound: 1.0
  starts: 6
  seed: 1
surface:
  r_values: [1, 1.25, 1.5, 2, 4, 8, 16, 32, 64]
  A_values: auto_around_mle
```

可选非主分析：

```yaml
species_cuts:
  topology_k: [16, 24]
```

### `species_tree.nwk`

- Newick、analysis-ready、rooted且严格二叉；
- tip名称为唯一`host_id`；
- tip集合等于included genes映射出的host集合；
- branch lengths允许存在，但v2统一按unit edge计算。

### `gene_tree.nwk`

- Newick，tip为唯一`gene_id`；
- tip集合等于`include=1`的gene集合；
- rootedness在YAML声明；
- 无根二叉树允许Newick根degree 3；
- 只用于成员、GeneCut和结构审计，不拆分D/T。

### `genes.tsv`

| 必需列 | 说明 |
|---|---|
| `gene_id` | gene-tree tip唯一ID |
| `host_id` | species-tree tip唯一ID |
| `include` | `1/true`进入分析 |

允许附加`contig_id,start,end,strand,protein_length,annotation_source,sequence_sha256,qc_status`等provenance列，但默认likelihood忽略。

## 2. 适用性状态

| 数据 | 状态 | 行为 |
|---|---|---|
| 没有C | `NOT_APPLICABLE_NO_COEXIST` | 不拟合free model |
| 没有S | `NOT_IDENTIFIABLE_NO_STATE_CONTRAST` | 不拟合r |
| S和C均有 | `FIT_ELIGIBLE` | 完整分析 |

## 3. 输出目录

### 输入与状态

| 文件 | 内容 |
|---|---|
| `input_audit.json` | 输入维度、SHA、edge policy |
| `validated_gene_to_host.tsv` | 验证后的mapping |
| `host_copy_counts.tsv` | host copy数和S/C状态 |
| `gene_tree_tip_audit.tsv` | gene tip mapping状态 |
| `species_tree_tip_audit.tsv` | host tip copy状态 |
| `applicability.json` | 是否适合拟合r |

### 正式检验与拟合

| 文件 | 内容 |
|---|---|
| `clustering_test.json` | Fitch与置换p值 |
| `coexistence_test.tsv` | 单行主结果 |
| `fit.json` | null/free、LRT、profile、事件期望 |
| `r_profile.tsv` | 固定r后profile优化A |
| `likelihood_surface.tsv` | `r x A` likelihood |
| `event_expectations_surface.tsv` | surface点宏观事件期望 |
| `branch_state_posteriors.tsv` | 每条边四状态posterior |
| `expected_macro_edge_transitions.tsv` | free MLE事件汇总 |

### Cut与汇总

| 文件 | 内容 |
|---|---|
| `gene_cut_audit.json` | gene split审计 |
| `gene_cut_manifest.tsv` | 非平凡gene-tree splits |
| `gene_block_members.tsv` | split较小侧成员 |
| `species_cut_manifest.tsv` | full及可选topology cuts |
| `host_block_members.tsv` | host-to-block mapping |
| `cut_sensitivity.tsv` | full及可选abstract fits |
| `summary.json` | 机器可读总汇 |
| `SUMMARY_COMPLETE.json` | 输入和summary SHA完成marker |
| `README.md` | 人工结果解释，不由程序自动生成 |

Surface表以`family_id,species_cut_id,r,A`为联合键；branch表以`parent,child`标识边。

## 4. 禁止的隐式输入处理

核心程序不会根据taxonomy猜测host、下载GTDB、从reference tree诱导subtree、自动裁剪gene-tree context、合并重复gene ID、静默删除未映射gene，或将零copy host并入single。这些操作必须在输入冻结前显式完成并保存provenance。
