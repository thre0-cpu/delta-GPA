#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
将平衡常数转换为标准吉布斯自由能变化并绘制相关性图
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

# 设置matplotlib支持中文字体
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

def convert_to_dg_values(k_values, temperature):
    """将平衡常数转换为标准吉布斯自由能变化 ΔG'° = -RT ln(K')"""
    R = 8.314462618e-3  # 气体常数，kJ/(mol*K)
    dg_values = -R * temperature * np.log(k_values)
    return dg_values

def plot_dg_correlation_and_calculate_errors():
    # 读取过滤后的数据
    df = pd.read_csv('TECRDB_test_filtered_full.csv')
    
    print(f'数据集总行数: {len(df)}')
    
    # 过滤掉实验值（K_prime）缺失的数据点，只保留有实验值和计算值的数据
    df_filtered = df.dropna(subset=['K_prime', 'calculated_K_prime'])
    
    print(f'有实验值和计算值的数据点数量: {len(df_filtered)}')
    
    # 提取数值
    experimental_k_values = df_filtered['K_prime'].values
    calculated_k_values = df_filtered['calculated_K_prime'].values
    temperatures = df_filtered['temperature'].values
    
    # 将平衡常数转换为吉布斯自由能变化
    experimental_dg_values = convert_to_dg_values(experimental_k_values, temperatures)
    calculated_dg_values = convert_to_dg_values(calculated_k_values, temperatures)
    
    # 计算各种误差指标
    mae = np.mean(np.abs(experimental_dg_values - calculated_dg_values))  # 平均绝对误差
    rmse = np.sqrt(np.mean((experimental_dg_values - calculated_dg_values)**2))  # 均方根误差
    mre = np.mean(np.abs((experimental_dg_values - calculated_dg_values) / (np.abs(experimental_dg_values) + 1e-10))) * 100  # 平均相对误差 (%)
    
    # 计算相关系数
    correlation, p_value = stats.pearsonr(experimental_dg_values, calculated_dg_values)
    
    print(f'吉布斯自由能变化 - 平均绝对误差 (MAE): {mae:.4f} kJ/mol')
    print(f'吉布斯自由能变化 - 均方根误差 (RMSE): {rmse:.4f} kJ/mol')
    print(f'吉布斯自由能变化 - 平均相对误差 (MRE): {mre:.2f}%')
    print(f'吉布斯自由能变化 - 皮尔逊相关系数: {correlation:.4f}')
    print(f'吉布斯自由能变化 - 相关性p值: {p_value:.2e}')
    
    # 计算决定系数R^2
    ss_res = np.sum((experimental_dg_values - calculated_dg_values) ** 2)
    ss_tot = np.sum((experimental_dg_values - np.mean(experimental_dg_values)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    print(f'吉布斯自由能变化 - 决定系数 R^2: {r_squared:.4f}')
    
    # 绘制吉布斯自由能变化的散点图
    plt.figure(figsize=(10, 8))
    
    # 绘制散点图
    plt.scatter(experimental_dg_values, calculated_dg_values, alpha=0.6, s=30)
    
    # 绘制y=x线（理想预测线）
    min_val = min(min(experimental_dg_values), min(calculated_dg_values))
    max_val = max(max(experimental_dg_values), max(calculated_dg_values))
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', label='y=x (理想预测)', linewidth=1)
    
    # 添加图表标题和标签
    plt.title('实验吉布斯自由能变化 vs 计算吉布斯自由能变化', fontsize=14)
    plt.xlabel('实验吉布斯自由能变化 ΔG\'°$_{exp}$ (kJ/mol)', fontsize=12)
    plt.ylabel('计算吉布斯自由能变化 ΔG\'°$_{calc}$ (kJ/mol)', fontsize=12)
    
    # 添加误差信息到图例
    plt.legend([f'y=x (理想预测)', 
                f'数据点 (n={len(df_filtered)})\nMAE={mae:.2f} kJ/mol\nRMSE={rmse:.2f} kJ/mol\nCorr={correlation:.3f}'],
               loc='upper left')
    
    # 添加网格
    plt.grid(True, alpha=0.3)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    plt.savefig('experimental_vs_calculated_dg_correlation.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 绘制误差直方图
    errors = experimental_dg_values - calculated_dg_values
    
    plt.figure(figsize=(10, 6))
    plt.hist(errors, bins=30, alpha=0.7, edgecolor='black')
    plt.xlabel('误差 ΔG\'°$_{exp}$ - ΔG\'°$_{calc}$ (kJ/mol)')
    plt.ylabel('频数')
    plt.title('ΔG\'° 误差分布')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('dg_error_distribution.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 计算一些额外的统计量
    print(f'\n额外统计信息:')
    print(f'实验ΔG\'°范围: {experimental_dg_values.min():.4f} - {experimental_dg_values.max():.4f} kJ/mol')
    print(f'计算ΔG\'°范围: {calculated_dg_values.min():.4f} - {calculated_dg_values.max():.4f} kJ/mol')
    print(f'实验ΔG\'°中位数: {np.median(experimental_dg_values):.4f} kJ/mol')
    print(f'计算ΔG\'°中位数: {np.median(calculated_dg_values):.4f} kJ/mol')
    print(f'实验ΔG\'°平均值: {np.mean(experimental_dg_values):.4f} kJ/mol')
    print(f'计算ΔG\'°平均值: {np.mean(calculated_dg_values):.4f} kJ/mol')

def main():
    plot_dg_correlation_and_calculate_errors()

if __name__ == '__main__':
    main()