# ATP + 葡萄糖反应热力学分析 - 技术报告

## 项目概述

本项目使用CC方法（组分贡献法）分析了ATP + 葡萄糖反应（'ATP + 葡萄糖 = 葡萄糖-6-磷酸 + ADP'）的热力学平衡特性，考察了不同pH和pMg条件对反应的影响。

## 核心技术方法

### 1. 组分贡献法 (Component Contribution, CC)
- 使用 `equilibrator_api` 库的 ComponentContribution 类
- 利用KEGG ID (ATP: C00002, 葡萄糖: C00031, 葡萄糖-6-磷酸: C00092, ADP: C00008) 进行化合物识别
- 计算在特定条件下的标准反应自由能变化 (ΔG'°)

### 2. 平衡计算方法
- 使用解析解方法求解平衡浓度
- 对于反应 ATP + 葡萄糖 = 葡萄糖-6-磷酸 + ADP
- 解析式: x = sqrt(K_eq) * [initial] / (1 + sqrt(K_eq))

### 3. 条件变量分析
- 温度: 298.15 K
- 离子强度: 0.25 M
- pH范围: 5.0 - 10.0
- pMg范围: 2.0 - 10.0
- 特定条件: pH 9, pMg 3.0

## 重要发现

### 1. 基本反应特性
- 在pH 9, pMg 3.0, 离子强度 0.25 M条件下:
  - ΔG'° = -28.85 kJ/mol
  - 平衡常数 K_eq = 1.13×10^5
  - ATP平衡浓度: 0.000003 M (几乎完全反应)

### 2. pMg影响分析
- ATP平衡浓度随pMg增加而降低
- 最低ATP浓度: 0.00000157 M (pMg=14.0)
- 最高ATP浓度: 0.00000401 M (pMg=2.0)
- 生理pMg(3.0)下ATP浓度: 0.00000296 M

### 3. 二维参数影响 (pH-pMg)
- 99个数据点的二维分析 (11个pH值 × 9个pMg值)
- ATP平衡浓度范围: 4.56×10⁻⁷ M - 8.88×10⁻⁵ M
- ΔG'°范围: -38.14 kJ/mol - -11.54 kJ/mol
- 生理条件估算(pH≈7.5, pMg≈3.0): [ATP]≈1.61×10⁻⁵ M, ΔG'°≈-20.40 kJ/mol

## 数据处理与可视化

### 1. 数据存储
- 生成了多个CSV数据文件
- 包括一维分析和二维参数分析的数据

### 2. 可视化方法
- 线图: ATP浓度 vs pMg
- 热图: ATP浓度和ΔG'°随pH和pMg变化的二维热图
- 使用matplotlib进行可视化，并解决了中文字体问题

## 关键代码模块

### 1. CC方法实现
```python
from equilibrator_api import ComponentContribution, Q_
cc = ComponentContribution()
cc.p_h = Q_(pH)
cc.p_mg = Q_(pMg)
cc.ionic_strength = Q_(f'{I}M')
cc.temperature = Q_(f'{T}K')
dg_prime = cc.standard_dg_prime(parsed_rxn)
```

### 2. 平衡计算
```python
K_eq = math.exp(-delta_g_prime_j_per_mol / (R * T))
x_solution = sqrt(K_eq) * 0.001 / (1 + sqrt(K_eq))
```

## 技术挑战与解决方案

1. **中文字体问题**: 使用 `plt.rcParams['font.sans-serif']` 设置中文字体
2. **f-string语法错误**: 通过变量替换解决反斜杠问题
3. **数据维度处理**: 正确重塑数据以适应热图可视化需求

## 生物学意义

- 揭示了pH和pMg对ATP相关反应的重要影响
- 为细胞内代谢调控提供了热力学基础
- 说明了细胞可通过调节pH和Mg²⁺浓度来精确控制代谢反应的平衡位置

## 项目文件结构

- `main.ipynb`: 核心工作流程和可视化
- `atp_concentration_vs_pMg.csv`: pMg影响数据
- `atp_concentration_vs_ph_pmg.csv`: 二维参数影响数据
- `atp_concentration_heatmap.png`: ATP浓度热图
- `delta_g_heatmap.png`: ΔG'°热图
- `temp_files/`: 存放临时脚本
- `TECH_report.md`: 本技术报告
- `CODE_report.md`: 代码分析报告
- `README.md`: 工作区说明