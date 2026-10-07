# coexist-gl 安装、部署、使用与结果读取

## 1. 环境

- Python >=3.10；
- Biopython >=1.80；
- NumPy >=1.24；
- SciPy >=1.10；
- PyYAML >=6.0。

这是Python package，不需要C/C++编译器。

## 2. 安装

Windows：

```powershell
cd <repository-root>
python -m pip install -e .
coexist-gl --help
```

Linux/macOS：

```bash
cd /path/to/coexistence_gain_loss_toolkit
python3 -m pip install -e .
coexist-gl --help
```

不安装运行：

```powershell
$env:PYTHONPATH = "$PWD\src"
python -m coexist_gl.cli --help
```

## 3. 测试

```powershell
$env:PYTHONPATH = "$PWD\src"
python -m unittest discover -s tests -v
```

预期5个测试全部`OK`。

## 4. 验证输入

```powershell
coexist-gl validate --config examples\actin_current584\input\family.yaml
```

成功时返回维度、single/coexist数、输入SHA和`status=PASS`，不进行拟合。

## 5. 完整运行

```powershell
coexist-gl run --config examples\actin_current584\input\family.yaml --output examples\actin_current584\results
```

actin示例通常为十余秒。结束后应存在`results/SUMMARY_COMPLETE.json`，且其中输入SHA应与当前input一致。

## 6. 推荐读结果顺序

1. `input_audit.json`：输入闭合；
2. `applicability.json`：family是否适用；
3. `host_copy_counts.tsv`：观测状态；
4. `clustering_test.json`：树上聚集；
5. `coexistence_test.tsv`：主结果；
6. `fit.json`：参数和事件期望；
7. `r_profile.tsv`：r区间是否闭合；
8. `branch_state_posteriors.tsv`：候选状态变化边；
9. surfaces：参数敏感性；
10. `SUMMARY_COMPLETE.json`：结果与输入匹配。

## 7. 批量families

每个family使用独立`input/`和`results/`目录。family之间可以并行，但不能共享output目录。合并`coexistence_test.tsv`后，应分别对`p_cluster`和`p_persistence_one_sided`做多重检验校正；程序目前不自动执行FDR。

## 8. 常见错误

- `gene tip/mapping mismatch`：gene tree与included genes不完全相同；
- `species tip/mapped-host mismatch`：species tree含零copy host或mapping含tree外host；
- `tree is not binary`：输入树未满足二叉合同；
- `R_NOT_IDENTIFIABLE`：r profile上界开放，不能把扫描最大值当置信上界。

## 9. 安全结果表述

可写：

> Coexistence states were phylogenetically clustered and a unit-edge state model allowing enhanced `C -> C` persistence improved the likelihood relative to the corresponding `r=1` collapsed gain-resolution model.

不可写本程序“识别了具体D/T/L次数”，因为输出是宏观S/C边转移期望。
