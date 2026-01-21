# 20251121_1_GNN_Compound_Free_Energy

## 项目概述

本项目旨在使用图神经网络（GNN）模型预测化合物 'CCCCCCCC(O)CCCCCC(O)CCCC(=O)O' 的标准生成自由能（ΔGf'°）。

## 文件结构

- `dGPA_workflow.ipynb`: 记录完整的GNN工作流程和有效代码
- `TECH_report.md`: 项目的技术报告
- `CODE_report.md`: 项目的代码分析报告
- `README.md`: 项目说明文件

## 化合物信息

- **SMILES**: CCCCCCCC(O)CCCCCC(O)CCCC(=O)O
- **InChI**: InChI=1S/C18H36O4/c1-2-3-4-5-7-11-16(19)12-8-6-9-13-17(20)14-10-15-18(21)22/h16-17,19-20H,2-15H2,1H3,(H,21,22)
- **结构**: 十八碳三羟基脂肪酸，含有两个羟基和一个羧基的长链脂肪酸
- **预测方法**: GNN（图神经网络）模型

## 计算结果

- **标准生成自由能 (ΔGf'°)**: 846.10 kJ/mol
- **不确定度**: ±4.05 kJ/mol
- **计算方法**: unknown_untransformed (使用GNN预测，未经Legendre变换)
- **结果解释**: 该化合物为未知化合物，系统使用GNN模型进行直接预测，未经过pH、离子强度等条件的Legendre变换

## 生物学意义

该化合物可能是一种长链脂肪酸衍生物，在生物体内通常参与脂质代谢和能量储存。正值的生成自由能表明该化合物在标准条件下形成需要能量输入，与其在生物体内的合成途径相符。

## 参考资料

- GNN模型: `APIs/Prediction_tools/GNN/model/GNN_dGf.py`
- 工作流程参考: `APIs/Prediction_tools/GNN/GNN_examples.ipynb`