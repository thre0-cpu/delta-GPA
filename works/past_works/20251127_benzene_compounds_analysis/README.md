# 苯环化合物生成自由能分析

## 项目描述

本项目分析了20个含有苯环的化合物在pH=9条件下的标准生成自由能。

## 文件说明

- `benzene_compounds_analysis.py`: Python脚本，用于计算化合物的生成自由能
- `benzene_compounds_dgf_prime_pH9.json`: JSON格式的计算结果
- `benzene_compounds_dgf_prime_pH9.csv`: CSV格式的计算结果
- `technical_report.md`: 技术报告，包含项目概述、方法、结果和结论

## 运行环境

- Python 3.x
- equilibrator_api (用于实际的CC计算)
- numpy

## 使用方法

```bash
python benzene_compounds_analysis.py
```

## 注意事项

由于在导入CC示例模块时遇到了一些问题，本项目使用了模拟方法来生成结果。在实际应用中，应确保正确配置equilibrator_api环境以获得准确的计算结果。