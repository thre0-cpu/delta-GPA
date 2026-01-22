# dGPredictor 工具包（Moiety Contribution 方法）

基于**分子签名（Molecular Signature）**和**贝叶斯岭回归**的反应自由能预测方法。

## 📁 文件结构

```
├── MC_examples.ipynb              # 使用示例和函数封装
├── README.md                      # 本文档
└── dGPredictor/                   # 核心模块
    ├── data/                      # 数据文件
    │   ├── cache_compounds_20160818.csv          # KEGG化合物SMILES数据库
    │   ├── decompose_vector_ac.json              # 分子签名 (radius=1)
    │   ├── decompose_vector_ac_r2_*.json         # 分子签名 (radius=2)
    │   ├── group_names_r1.txt                    # 基团名称列表 (r=1)
    │   ├── group_names_r2_*.txt                  # 基团名称列表 (r=2)
    │   └── Test_KEGG_all_grp.mat                 # 训练数据
    ├── model/
    │   └── M12_model_BR.pkl       # 预训练贝叶斯岭回归模型
    └── generate_model.py          # 模型生成脚本
```

## 📄 核心原理

1. **分子签名分解**：将分子按 ECFP-like 方式分解为子结构特征向量（radius=1 和 radius=2）
2. **反应规则生成**：计算反应前后分子签名的变化
3. **贝叶斯岭回归**：使用预训练模型预测 ΔrG'°
4. **不确定度估计**：贝叶斯后验分布提供置信区间

## 🚀 快速开始

```python
# 在 MC_examples.ipynb 中已封装好函数

# 使用 KEGG ID 预测
dG, std = predict_dG("C00002 + C00001 <=> C00008 + C00009")  # ATP水解
print(f"ΔrG'° = {dG:.2f} ± {std:.2f} kJ/mol")

# 使用 SMILES 预测（支持新分子）
dG, std = predict_dG_from_smiles(
    reactants={"CCO": 1},      # 乙醇
    products={"CC=O": 1, "O": 1}  # 乙醛 + 水
)
```

## 📊 与其他方法对比

| 特性 | CC (eQuilibrator) | MC (dGPredictor) | GNN |
|------|------------------|------------------|-----|
| 分解方法 | 基团贡献 | 分子签名 | 图神经网络 |
| 需要 pKa | 是 (ChemAxon) | 否 | 否 |
| 新分子支持 | 需计算 pKa | 只需 SMILES | 只需 SMILES |
| 不确定度 | 协方差矩阵 | 贝叶斯后验 | 集成/MC Dropout |

## 📚 参考文献

- Du, J. et al. "dGPredictor: Automated fragmentation method for metabolic reaction free energy prediction"
- GitHub: https://github.com/maranasgroup/dGPredictor

## ⚠️ 注意事项

- `model/` 文件夹内为预训练模型，**不要修改**
- 如模型丢失，运行 `python dGPredictor/generate_model.py` 重新生成
