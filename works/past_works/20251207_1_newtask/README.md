# ATP + 葡萄糖反应热力学分析项目

## 项目简介

本项目使用CC方法（组分贡献法）分析了ATP + 葡萄糖反应（'ATP + 葡萄糖 = 葡萄糖-6-磷酸 + ADP'）的热力学平衡特性，并考察了不同pH和pMg条件对反应的影响。

## 项目目标

1. 计算在特定条件下（298.15 K, pH 9, pMg 3.0, 离子强度 0.25 M）反应的平衡组分
2. 分析pMg对ATP平衡浓度的影响
3. 扩展分析pH和pMg对反应的综合影响（二维参数空间）

## 项目结构

```
├── main.ipynb                 # 核心工作流程和可视化
├── temp_files/                # 临时脚本文件夹
│   ├── equilibrium_calculation.py
│   ├── pMg_analysis.py
│   ├── ph_pmg_analysis.py
│   ├── plot_atp_vs_pmg.py
│   ├── plot_atp_vs_pmg_chinese.py
│   ├── test_cc.py
│   ├── test_cc_v2.py
│   ├── comprehensive_report.md
│   ├── comprehensive_report_2d.md
│   ├── pMg_analysis.ipynb
│   ├── pMg_analysis_report.md
│   └── report.md
├── results/                   # 结果文件夹
│   ├── atp_concentration_vs_pMg.csv
│   ├── atp_concentration_vs_ph_pmg.csv
│   ├── atp_concentration_heatmap.png
│   ├── atp_concentration_vs_pMg.png
│   ├── atp_concentration_vs_pMg_plot.png
│   ├── atp_concentration_vs_pMg_plot_chinese_font.png
│   ├── delta_g_heatmap.png
│   └── ...
├── reports/                   # 报告文件夹
│   ├── TECH_report.md         # 技术报告
│   └── CODE_report.md         # 代码分析报告
└── README.md                  # 本说明文件
```

## 核心发现

### 基本反应特性
- 在pH 9, pMg 3.0, 离子强度 0.25 M条件下:
  - ΔG'° = -28.85 kJ/mol
  - 平衡常数 K_eq = 1.13×10^5
  - ATP平衡浓度: 0.000003 M (几乎完全反应)

### pMg影响
- ATP平衡浓度随pMg增加而降低
- 生理pMg(3.0)下ATP浓度: 0.00000296 M

### 二维参数影响 (pH-pMg)
- 99个数据点的二维分析 (pH 5.0-10.0, pMg 2.0-10.0)
- ATP平衡浓度范围: 4.56×10⁻⁷ M - 8.88×10⁻⁵ M
- 生理条件估算(pH≈7.5, pMg≈3.0): [ATP]≈1.61×10⁻⁵ M

## 使用方法

### 环境要求
- Python 3.8+
- equilibrator_api
- numpy, pandas, matplotlib
- scipy

### 运行流程
1. `main.ipynb` 包含了完整的分析流程，可在Jupyter环境中运行
2. 所有图表和数据均已生成，可直接查看结果

## 生物学意义

本项目揭示了pH和pMg对ATP + 葡萄糖反应的重要影响，为理解细胞内代谢调控机制提供了热力学基础。结果表明细胞可通过调节pH和Mg²⁺浓度来精确控制代谢反应的平衡位置。

## 技术要点

- 使用equilibrator_api的ComponentContribution类进行热力学计算
- 通过KEGG ID识别化合物 (ATP: C00002, 葡萄糖: C00031, 葡萄糖-6-磷酸: C00092, ADP: C00008)
- 解析求解平衡浓度，避免数值计算误差
- 生成热图可视化二维参数空间的影响

## 作者与时间

- 项目创建时间: 2025年12月7日
- 目录: 20251207_1_New_Task