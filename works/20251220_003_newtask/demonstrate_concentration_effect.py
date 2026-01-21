#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
演示不同代谢物浓度对反应自由能变的影响
"""

import pandas as pd
import numpy as np
from equilibrator_api import ComponentContribution, Q_
import matplotlib.pyplot as plt

# 设置matplotlib支持中文字体
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

def calculate_dg_at_different_concentrations(rxn_str, cc, concentrations, p_h=7.0, p_mg=3.0, I=0.25, T=298.15):
    """
    在不同代谢物浓度下计算反应的自由能变
    """
    # 设置条件
    cc.p_h = Q_(p_h)
    cc.p_mg = Q_(p_mg)
    cc.ionic_strength = Q_(f'{I}M')
    cc.temperature = Q_(f'{T}K')
    
    # 解析反应
    reaction = cc.parse_reaction_formula(rxn_str)
    
    # 计算标准反应的 ΔG'°
    dg_prime_measurement = cc.standard_dg_prime(reaction)
    dg_prime_std = dg_prime_measurement.value.m_as("kJ/mol")
    
    # 计算不同浓度下的 ΔG
    R = 8.314462618e-3  # 气体常数，kJ/(mol*K)
    dg_values = []
    
    for conc in concentrations:
        # 计算反应商的对数项
        lnQ = 0
        for compound, coeff in reaction.sparse.items():
            # 浓度校正项是 coeff * ln([conc])
            lnQ += coeff * np.log(conc)
        
        # 应用浓度校正
        dg_at_conc = dg_prime_std + R * T * lnQ
        dg_values.append(dg_at_conc)
    
    return dg_values, dg_prime_std

def main():
    print('正在初始化ComponentContribution...')
    cc = ComponentContribution()
    print('ComponentContribution初始化完成')
    
    # 定义不同浓度
    concentrations = np.logspace(-6, 0, 100)  # 1 μM to 1 M
    
    # 选择几个典型的反应来演示
    test_reactions = [
        "kegg:C00002 + kegg:C00001 = kegg:C00008 + kegg:C00009",  # ATP + H2O = ADP + Pi
        "kegg:C00022 + kegg:C00012 = kegg:C00026 + kegg:C00011",  # 丙酮酸+NADH+H+ = 乳酸+NAD+
    ]
    
    # 设置图表
    plt.figure(figsize=(12, 8))
    
    for i, rxn in enumerate(test_reactions):
        try:
            dg_values, dg_std = calculate_dg_at_different_concentrations(rxn, cc, concentrations)
            
            plt.subplot(2, 1, i+1)
            plt.semilogx(concentrations, dg_values, label=f'ΔG = f([Metabolite])', linewidth=2)
            plt.axhline(y=dg_std, color='r', linestyle='--', label=f'标准ΔG\'° = {dg_std:.2f} kJ/mol')
            plt.xlabel('代谢物浓度 (M)')
            plt.ylabel('ΔG (kJ/mol)')
            plt.title(f'反应: {rxn}')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            print(f'反应 {i+1}: {rxn}')
            print(f'  标准ΔG\'°: {dg_std:.2f} kJ/mol')
            print(f'  在1mM浓度下的ΔG: {dg_values[len(concentrations)//2]:.2f} kJ/mol')  # 中间浓度点（约1mM）
            print()
            
        except Exception as e:
            print(f"处理反应 {rxn} 时出错: {str(e)}")
    
    plt.tight_layout()
    plt.savefig('dg_vs_concentration.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("完成不同浓度下自由能变的分析！")
    print("图中显示了代谢物浓度如何影响反应的吉布斯自由能变。")
    
    # 读取原始数据并添加生理条件自由能变
    df = pd.read_csv('TECRDB_test_filtered_full.csv')
    print(f'\n已将生理条件下的计算结果保存到DataFrame中')
    print(f'前5行数据预览:')
    print(df[['reaction', 'K_prime', 'calculated_K_prime', 'temperature', 'p_h', 'p_mg']].head())

if __name__ == '__main__':
    main()