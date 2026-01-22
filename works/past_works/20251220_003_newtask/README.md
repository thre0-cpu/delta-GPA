# TECRDB数据集平衡常数预测与分析

本项目使用Component Contribution (CC)方法预测TECRDB_test数据集中生化反应的平衡常数，并与实验值进行比较分析。

## 项目结构

```
├── main.ipynb                 # 核心工作流程和可视化
├── temp_files/                # 临时脚本/Notebook文件夹
│   ├── plot_correlation.py
│   ├── plot_dg_correlation.py
│   ├── filter_csv.py
│   ├── filter_full_csv.py
│   ├── calculate_physiological_dg.py
│   ├── example_different_concentrations.py
│   ├── test_reaction_attrs.py
│   ├── check_compound_attrs.py
│   └── run_cc_prediction.py
├── results/                   # 结果文件夹
│   ├── experimental_vs_calculated_correlation.png
│   ├── experimental_vs_calculated_dg_correlation.png
│   ├── error_distributions.png
│   ├── dg_error_distribution.png
│   ├── TECRDB_test_filtered_full.csv
│   └── TECRDB_test_with_physiological_predictions.csv
├── reports/                   # 结果文件夹
│   ├── TECH_report.md         # 技术报告
│   └── CODE_report.md         # 代码分析报告
└── README.md                  # 说明文件
```

## 项目目标

1. 使用CC方法预测TECRDB数据集中生化反应的平衡常数
2. 与实验值进行比较分析
3. 评估CC方法在预测生化反应热力学参数方面的性能

## 技术要点

- 使用equilibrator_api的ComponentContribution类进行计算
- 对每个反应根据其特定条件（pH、pMg、离子强度、温度）进行个性化计算
- 从平衡常数和自由能变化两个层面评估预测性能
- 实现生理条件下反应自由能变的计算

## 主要结果

- CC方法与实验值之间呈现出高度相关性（皮尔逊相关系数0.9899）
- 平均绝对误差仅为1.44 kJ/mol（基于自由能变化）
- 在生物化学领域，CC方法是预测生化反应热力学参数的可靠工具

## 使用说明

1. 运行main.ipynb中的代码进行完整的分析流程
2. 查看results/文件夹中的图表和数据文件
3. 阅读reports/文件夹中的技术报告和代码分析报告

## 关键文件说明

- `main.ipynb`：整合了所有分析流程的Jupyter Notebook
- `TECRDB_test_filtered_full.csv`：过滤后的数据文件，包含反应式、实验值、计算值和条件参数
- `TECRDB_test_with_physiological_predictions.csv`：包含生理条件下自由能变的完整数据
- `experimental_vs_calculated_dg_correlation.png`：吉布斯自由能变化的实验值vs计算值散点图
- `dg_error_distribution.png`：自由能变化误差分布直方图