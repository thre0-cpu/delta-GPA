# 不同pMg下ATP + 葡萄糖反应中ATP浓度分析
import subprocess
import json
import sys
import os
from typing import Optional, Tuple, Union

import numpy as np
import numpy.typing as npt
from rdkit import Chem

from equilibrator_api import ComponentContribution, Q_
from scipy.optimize import fsolve
from scipy.constants import R
import math
import matplotlib
matplotlib.use('Agg')  # 使用非GUI后端
import matplotlib.pyplot as plt
import pandas as pd

print("正在初始化ComponentContribution实例...")

# 初始化CC类
cc = ComponentContribution()
print("ComponentContribution实例已创建")

# 定义pMg范围
pMg_values = np.arange(1.0, 14.1, 0.5)
atp_concentrations = []

# 固定条件
p_h = 9.0
I = 0.25  # 离子强度 (M)
T = 298.15  # 温度 (K)

# 定义反应，使用KEGG ID
reaction_formula = 'C00002 + C00031 = C00092 + C00008'
parsed_rxn = cc.parse_reaction_formula(reaction_formula)
print('反应解析完成')

# 遍历不同的pMg值
for i, pMg in enumerate(pMg_values):
    # 设置反应条件
    cc.p_h = Q_(p_h)
    cc.p_mg = Q_(pMg)
    cc.ionic_strength = Q_(f'{I}M')
    cc.temperature = Q_(f'{T}K')
    
    # 计算标准反应自由能变化
    dg_prime = cc.standard_dg_prime(parsed_rxn)
    
    # 计算平衡常数 K_eq
    delta_g_prime_kj_per_mol = dg_prime.value.m  # kJ/mol
    delta_g_prime_j_per_mol = delta_g_prime_kj_per_mol * 1000  # J/mol
    RT = R * T
    ln_keq = -delta_g_prime_j_per_mol / RT
    K_eq = math.exp(ln_keq)
    
    # 计算平衡时ATP浓度
    # 解析解: sqrt(K_eq) = x / (0.001-x)
    # x = sqrt(K_eq) * 0.001 / (1 + sqrt(K_eq))
    sqrt_keq = math.sqrt(K_eq)
    x_solution = sqrt_keq * 0.001 / (1 + sqrt_keq)
    
    # 平衡时ATP浓度
    atp_eq = 0.001 - x_solution
    
    atp_concentrations.append(atp_eq)
    
    # 输出进度（每5个值输出一次）
    if i % 5 == 0:
        print(f'pMg: {pMg:.1f}, ATP浓度: {atp_eq:.8f} M, K_eq: {K_eq:.2e}')

print(f'完成计算，共处理 {len(pMg_values)} 个pMg值')

# 绘制ATP浓度随pMg变化的图表
plt.figure(figsize=(12, 8))
plt.plot(pMg_values, atp_concentrations, 'b-o', markersize=4, linewidth=1.5)
plt.xlabel('pMg')
plt.ylabel('ATP平衡浓度 (M)')
plt.title('不同pMg下ATP的平衡浓度 (298.15 K, pH 9, 0.25 M 离子强度)')
plt.grid(True, alpha=0.3)
plt.yscale('log')  # 使用对数刻度更好地显示浓度变化

# 在对数刻度上添加文本标注关键点
min_idx = atp_concentrations.index(min(atp_concentrations))
max_idx = atp_concentrations.index(max(atp_concentrations))

plt.annotate(f'Min: ({pMg_values[min_idx]:.1f}, {min(atp_concentrations):.2e})',
             xy=(pMg_values[min_idx], min(atp_concentrations)),
             xytext=(pMg_values[min_idx]+1.5, min(atp_concentrations)*2),
             arrowprops=dict(arrowstyle='->', color='red', lw=0.8))

plt.annotate(f'Max: ({pMg_values[max_idx]:.1f}, {max(atp_concentrations):.2e})',
             xy=(pMg_values[max_idx], max(atp_concentrations)),
             xytext=(pMg_values[max_idx]-2, max(atp_concentrations)/2),
             arrowprops=dict(arrowstyle='->', color='red', lw=0.8))

# 保存图表
plt.savefig('atp_concentration_vs_pMg.png', dpi=300, bbox_inches='tight')
print('图表已保存为atp_concentration_vs_pMg.png')

# 创建详细数据表
data = {
    'pMg': pMg_values,
    'ATP平衡浓度 (M)': atp_concentrations
}

df = pd.DataFrame(data)
print('\n关键数据点:')
print(df)

# 保存数据到CSV文件
df.to_csv('atp_concentration_vs_pMg.csv', index=False)
print('\n数据已保存到atp_concentration_vs_pMg.csv')

# 打印总结信息
print(f'\n总结:')
print(f'pMg=3 (生理条件) 时，ATP浓度: {atp_concentrations[int((3.0-1.0)/0.5)]:.8f} M')
print(f'最低ATP浓度: {min(atp_concentrations):.8f} M (在pMg={pMg_values[atp_concentrations.index(min(atp_concentrations))]}处)')
print(f'最高ATP浓度: {max(atp_concentrations):.8f} M (在pMg={pMg_values[atp_concentrations.index(max(atp_concentrations))]}处)')
print(f'ATP浓度变化范围: {min(atp_concentrations):.8f} M - {max(atp_concentrations):.8f} M')