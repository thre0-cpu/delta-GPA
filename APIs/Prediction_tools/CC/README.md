# eQuilibrator 工具包

这个文件夹包含了使用 eQuilibrator API 进行生化反应热力学分析的工具和文档。

## 📁 文件结构

```
├── CC_examples.ipynb          # Component Contribution 方法示例代码
├── eQuilibrator_API.md        # eQuilibrator API 官方参考文档
└── eQuilibrator_references.md # eQuilibrator 中文使用指南
```

## 📄 文件说明

### `CC_examples.ipynb`
- 包含基于组分贡献法（Component Contribution, CC）的函数示例
- 提供生成自由能预测（Δ_fG）和反应自由能预测（Δ_rG）的实用函数
- 可直接调用进行简单的生化反应热力学分析

### `eQuilibrator_API.md`
- eQuilibrator API 官方参考文档
- 详细的模块、类和函数说明
- 包含完整的 API 接口定义

### `eQuilibrator_references.md`
- 由 Claude 4.5 Sonnet 编写的中文参考指南
- 包含快速开始教程、核心模块介绍
- 提供完整的使用示例和最佳实践
- 涵盖错误处理和性能优化建议

## 🚀 快速开始

```python
from equilibrator_api import ComponentContribution

# 初始化
cc = ComponentContribution()

# 搜索化合物并计算反应自由能
rxn = cc.parse_reaction_formula("glucose + ATP => glucose-6-phosphate + ADP")
dg_prime = cc.standard_dg_prime(rxn)
print(f"ΔG'° = {dg_prime}")
```

## 📚 推荐阅读顺序

1. **新手入门**：先阅读 `eQuilibrator_references.md` 的"快速开始"部分
2. **实践操作**：参考 `CC_examples.ipynb` 中的函数示例
3. **深入学习**：查阅 `eQuilibrator_API.md` 了解完整 API
