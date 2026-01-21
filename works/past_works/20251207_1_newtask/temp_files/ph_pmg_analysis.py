# 同时考察pH和pMg对ATP + 葡萄糖反应的影响
import subprocess
import json
import sys
import os
from typing import Optional, Tuple, Union

import numpy as np
import numpy.typing as npt
from rdkit import Chem

from equilibrator_api import ComponentContribution, Q_
from scipy.constants import R
import math
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import font_manager

# 解决中文字体问题
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 初始化CC类
print("正在初始化ComponentContribution实例...")
cc = ComponentContribution()
print("ComponentContribution实例已创建")

# 定义pH和pMg范围
pH_values = np.arange(5.0, 10.1, 0.5)  # pH 5.0 到 10.0
pMg_values = np.arange(2.0, 10.1, 1.0)  # pMg 2.0 到 10.0

# 固定其他条件
I = 0.25  # 离子强度 (M)
T = 298.15  # 温度 (K)

# 定义反应，使用KEGG ID
reaction_formula = 'C00002 + C00031 = C00092 + C00008'
parsed_rxn = cc.parse_reaction_formula(reaction_formula)
print('反应解析完成')

# 为所有pH和pMg组合计算ATP平衡浓度
atp_concentrations = []
data_points = []

for pH in pH_values:
    for pMg in pMg_values:
        # 设置反应条件
        cc.p_h = Q_(pH)
        cc.p_mg = Q_(pMg)
        cc.ionic_strength = Q_(f'{I}M')
        cc.temperature = Q_(f'{T}K')
        
        # 计算标准反应自由能变化
        try:
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
            data_points.append({
                'pH': pH,
                'pMg': pMg,
                'ATP平衡浓度 (M)': atp_eq,
                'K_eq': K_eq,
                'ΔG\'° (kJ/mol)': delta_g_prime_kj_per_mol
            })
            
        except Exception as e:
            print(f"在pH={pH}, pMg={pMg}时计算出错: {e}")
            atp_concentrations.append(float('nan'))
            data_points.append({
                'pH': pH,
                'pMg': pMg,
                'ATP平衡浓度 (M)': float('nan'),
                'K_eq': float('nan'),
                'ΔG\'° (kJ/mol)': float('nan')
            })

# 创建数据框
df = pd.DataFrame(data_points)
print(f"计算完成，共{len(df)}个数据点")

# 保存数据到CSV
df.to_csv('atp_concentration_vs_ph_pmg.csv', index=False, encoding='utf-8')
print("数据已保存到atp_concentration_vs_ph_pmg.csv")

# 创建pH-pMg二维热图
fig, ax = plt.subplots(figsize=(14, 10))

# 重塑数据以适应热图
Z = np.array(atp_concentrations).reshape(len(pH_values), len(pMg_values))

# 创建热图
im = ax.imshow(Z, cmap='viridis', aspect='auto', extent=[pMg_values.min(), pMg_values.max(), pH_values.min(), pH_values.max()], origin='lower', vmin=np.nanmin(Z), vmax=np.nanmax(Z))

# 添加颜色条
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('ATP平衡浓度 (M)', fontsize=12)

# 设置坐标轴
ax.set_xlabel('pMg', fontsize=12)
ax.set_ylabel('pH', fontsize=12)
ax.set_title('ATP平衡浓度随pH和pMg变化的热图\n(298.15 K, 0.25 M 离子强度)', fontsize=14)

# 在图中添加数值标签（仅在较小的数据集上建议这样做）
if len(pH_values) <= 10 and len(pMg_values) <= 10:  # 避免标签过多
    for i, pH in enumerate(pH_values):
        for j, pMg in enumerate(pMg_values):
            text = ax.text(pMg, pH, f'{Z[i, j]:.2e}',
                           ha="center", va="center", color="white", fontsize=8)

plt.tight_layout()
plt.savefig('atp_concentration_heatmap.png', dpi=300, bbox_inches='tight')
print("热图已保存为atp_concentration_heatmap.png")
plt.show()

# 创建ΔG'°随pH和pMg变化的热图
fig, ax = plt.subplots(figsize=(14, 10))

# 获取ΔG'°数据
dG_values = df['ΔG\'° (kJ/mol)'].values.reshape(len(pH_values), len(pMg_values))

# 创建ΔG'°热图
im2 = ax.imshow(dG_values, cmap='RdBu_r', aspect='auto', extent=[pMg_values.min(), pMg_values.max(), pH_values.min(), pH_values.max()], origin='lower')

# 添加颜色条
cbar2 = plt.colorbar(im2, ax=ax)
cbar2.set_label('ΔG\'° (kJ/mol)', fontsize=12)

# 设置坐标轴
ax.set_xlabel('pMg', fontsize=12)
ax.set_ylabel('pH', fontsize=12)
ax.set_title('标准反应自由能变化随pH和pMg变化的热图\n(298.15 K, 0.25 M 离子强度)', fontsize=14)

# 在图中添加数值标签（仅在较小的数据集上建议这样做）
if len(pH_values) <= 10 and len(pMg_values) <= 10:
    for i, pH in enumerate(pH_values):
        for j, pMg in enumerate(pMg_values):
            text = ax.text(pMg, pH, f'{dG_values[i, j]:.1f}',
                           ha="center", va="center", color="white", fontsize=8)

plt.tight_layout()
plt.savefig('delta_g_heatmap.png', dpi=300, bbox_inches='tight')
print("ΔG'°热图已保存为delta_g_heatmap.png")
plt.show()

# 定义列名变量以避免f-string中的反斜杠问题
dg_column = 'ΔG\'° (kJ/mol)'

# 打印统计信息
print(f"\nATP平衡浓度统计:")
print(f"最小值: {df['ATP平衡浓度 (M)'].min():.2e} M")
print(f"最大值: {df['ATP平衡浓度 (M)'].max():.2e} M")
print(f"平均值: {df['ATP平衡浓度 (M)'].mean():.2e} M")

print(f"\nΔG'°统计:")
print(f"最小值: {df[dg_column].min():.2f} kJ/mol")
print(f"最大值: {df[dg_column].max():.2f} kJ/mol")
print(f"平均值: {df[dg_column].mean():.2f} kJ/mol")

# 检查在生理条件下的值（pH=7.4, pMg=3.0）
physiological_pH_idx = (df['pH'] - 7.4).abs().idxmin()
physiological_pMg_idx = (df['pMg'] - 3.0).abs().idxmin()

# 找到最接近生理条件的组合
physio_row = df.iloc[(abs(df['pH'] - 7.4) + abs(df['pMg'] - 3.0)).idxmin()]
print(f"\n生理条件(pH≈{physio_row['pH']:.1f}, pMg≈{physio_row['pMg']:.1f})下的ATP平衡浓度: {physio_row['ATP平衡浓度 (M)']:.2e} M")
print(f"生理条件下的ΔG'°: {physio_row[dg_column]:.2f} kJ/mol")