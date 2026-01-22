#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
绘制实验值与计算值的散点图并计算误差
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

# 设置matplotlib支持中文字体
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

def plot_correlation_and_calculate_errors():
    # 读取过滤后的数据
    df = pd.read_csv('TECRDB_test_filtered_full.csv')
    
    print(f'数据集总行数: {len(df)}')
    
    # 过滤掉实验值（K_prime）缺失的数据点，只保留有实验值和计算值的数据
    df_filtered = df.dropna(subset=['K_prime', 'calculated_K_prime'])
    
    print(f'有实验值和计算值的数据点数量: {len(df_filtered)}')
    
    # 计算误差
    experimental_values = df_filtered['K_prime'].values
    calculated_values = df_filtered['calculated_K_prime'].values
    
    # 计算各种误差指标
    mae = np.mean(np.abs(experimental_values - calculated_values))  # 平均绝对误差
    rmse = np.sqrt(np.mean((experimental_values - calculated_values)**2))  # 均方根误差
    mre = np.mean(np.abs((experimental_values - calculated_values) / experimental_values)) * 100  # 平均相对误差 (%)
    
    # 计算相关系数
    correlation, p_value = stats.pearsonr(experimental_values, calculated_values)
    
    print(f'平均绝对误差 (MAE): {mae:.4f}')
    print(f'均方根误差 (RMSE): {rmse:.4f}')
    print(f'平均相对误差 (MRE): {mre:.2f}%')
    print(f'皮尔逊相关系数: {correlation:.4f}')
    print(f'相关性p值: {p_value:.2e}')
    
    # 绘制散点图
    plt.figure(figsize=(10, 8))
    
    # 绘制散点图
    plt.scatter(experimental_values, calculated_values, alpha=0.6, s=30)
    
    # 绘制y=x线（理想预测线）
    min_val = min(min(experimental_values), min(calculated_values))
    max_val = max(max(experimental_values), max(calculated_values))
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', label='y=x (理想预测)', linewidth=1)
    
    # 设置对数刻度（因为平衡常数范围可能很大）
    plt.xscale('log')
    plt.yscale('log')
    
    # 添加图表标题和标签
    plt.title('实验平衡常数 vs 计算平衡常数\n(对数尺度)', fontsize=14)
    plt.xlabel('实验平衡常数 (K\'$_{exp}$)', fontsize=12)
    plt.ylabel('计算平衡常数 (K\'$_{calc}$)', fontsize=12)
    
    # 添加误差信息到图例
    plt.legend([f'y=x (理想预测)', 
                f'数据点 (n={len(df_filtered)})\nMAE={mae:.4f}\nRMSE={rmse:.4f}\nCorr={correlation:.4f}'],
               loc='upper left')
    
    # 添加网格
    plt.grid(True, alpha=0.3)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    plt.savefig('experimental_vs_calculated_correlation.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 绘制误差直方图
    errors = experimental_values - calculated_values
    relative_errors = (experimental_values - calculated_values) / experimental_values * 100
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # 绝对误差分布
    ax1.hist(errors, bins=30, alpha=0.7, edgecolor='black')
    ax1.set_xlabel('绝对误差 (K\'$_{exp}$ - K\'$_{calc}$)')
    ax1.set_ylabel('频数')
    ax1.set_title('绝对误差分布')
    ax1.grid(True, alpha=0.3)
    
    # 相对误差分布
    ax2.hist(relative_errors, bins=30, alpha=0.7, edgecolor='black')
    ax2.set_xlabel('相对误差 (%)')
    ax2.set_ylabel('频数')
    ax2.set_title('相对误差分布')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('error_distributions.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 计算一些额外的统计量
    print(f'\n额外统计信息:')
    print(f'实验值范围: {experimental_values.min():.4f} - {experimental_values.max():.2e}')
    print(f'计算值范围: {calculated_values.min():.4f} - {calculated_values.max():.2e}')
    print(f'实验值中位数: {np.median(experimental_values):.4f}')
    print(f'计算值中位数: {np.median(calculated_values):.4f}')
    print(f'实验值平均值: {np.mean(experimental_values):.4f}')
    print(f'计算值平均值: {np.mean(calculated_values):.4f}')
    
    # 计算决定系数R^2
    ss_res = np.sum((experimental_values - calculated_values) ** 2)
    ss_tot = np.sum((experimental_values - np.mean(experimental_values)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    print(f'决定系数 R^2: {r_squared:.4f}')

def main():
    plot_correlation_and_calculate_errors()

if __name__ == '__main__':
    main()