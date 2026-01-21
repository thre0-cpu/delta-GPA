#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
演示如何为Reaction对象中的每个物质分别设置浓度
"""

import numpy as np
from equilibrator_api import ComponentContribution, Q_, Reaction
from scipy.constants import R

def example_reaction_with_different_concentrations():
    """
    演示如何为不同物质设置不同浓度来计算非标准条件下的反应自由能
    """
    print('正在初始化ComponentContribution...')
    cc = ComponentContribution()
    print('ComponentContribution初始化完成')
    
    # 解析一个简单的反应
    rxn_str = "kegg:C00002 + kegg:C00001 = kegg:C00008 + kegg:C00009"  # ATP + H2O = ADP + Pi
    reaction = cc.parse_reaction_formula(rxn_str)
    
    print(f'反应: {rxn_str}')
    print('反应物和产物的化学计量系数:')
    for compound, coeff in reaction.sparse.items():
        # 尝试获取化合物名称
        compound_name = "Unknown"
        for identifier in compound.identifiers:
            if identifier.registry.namespace in ['kegg', 'bigg.metabolite', 'chebi', 'hmdb']:
                compound_name = identifier.accession
                break
        print(f'  {compound_name}: {coeff}')
    
    # 设置标准条件参数
    pH = 7.0
    pMg = 3.0
    ionic_strength = Q_('0.25M')
    temperature = Q_('298.15K')
    
    # 设置CC对象的条件
    cc.p_h = Q_(pH)
    cc.p_mg = Q_(pMg)
    cc.ionic_strength = ionic_strength
    cc.temperature = temperature
    
    # 计算标准条件下的ΔG'°
    dg_prime_standard = cc.standard_dg_prime(reaction)
    dg_std_kjmol = dg_prime_standard.value.m_as('kJ/mol')
    
    print(f'\\n标准条件下的ΔG\'°: {dg_std_kjmol:.2f} kJ/mol')
    
    # 现在演示如何为每个物质设置不同的浓度
    # 创建一个浓度字典，键是compound对象，值是浓度（单位M）
    concentrations = {}  # 单位：M
    for compound, coeff in reaction.sparse.items():
        # 获取化合物名称
        compound_name = "Unknown"
        for identifier in compound.identifiers:
            if identifier.registry.namespace in ['kegg', 'bigg.metabolite', 'chebi', 'hmdb']:
                compound_name = identifier.accession
                break

        if 'ATP' in compound_name or compound_name == 'ATP':
            concentrations[compound] = 1e-3  # 1 mM ATP
        elif 'H2O' in compound_name or compound_name == 'H2O':
            concentrations[compound] = 55.5  # 水的浓度 ~55.5 M
        elif 'ADP' in compound_name or compound_name == 'ADP':
            concentrations[compound] = 0.1e-3  # 0.1 mM ADP
        elif 'phosphate' in compound_name.lower() or 'Pi' in compound_name:
            concentrations[compound] = 5e-3   # 5 mM 磷酸盐
        else:
            concentrations[compound] = 1e-3   # 默认1 mM
    
    print('\\n各物质的浓度设置:')
    for compound, conc in concentrations.items():
        # 获取化合物名称
        compound_name = "Unknown"
        for identifier in compound.identifiers:
            if identifier.registry.namespace in ['kegg', 'bigg.metabolite', 'chebi', 'hmdb']:
                compound_name = identifier.accession
                break
        print(f'  {compound_name}: {conc} M ({conc*1000:.3f} mM)')
    
    # 计算反应商 (Q)
    ln_reaction_quotient = 0
    for compound, coeff in reaction.sparse.items():
        conc = concentrations[compound]
        ln_reaction_quotient += coeff * np.log(conc)  # coeff是化学计量系数
    
    print(f'\\n反应商的对数 ln(Q): {ln_reaction_quotient:.4f}')
    
    # 计算非标准条件下的ΔG
    RT = R * 1e-3 * float(temperature.magnitude)  # R*T，转换为kJ/mol单位
    dg_non_standard = dg_std_kjmol + RT * ln_reaction_quotient
    
    print(f'\\n非标准条件下的ΔG: {dg_non_standard:.2f} kJ/mol')
    print(f'标准条件到非标准条件的校正: {RT * ln_reaction_quotient:.2f} kJ/mol')
    
    # 计算反应的平衡常数
    K_prime = np.exp(-dg_std_kjmol / RT)
    print(f'\\n平衡常数 K\' : {K_prime:.2e}')
    
    # 计算反应的驱动力（使用实际浓度）
    Q_actual = np.exp(ln_reaction_quotient)
    delta_g_actual = dg_std_kjmol + RT * ln_reaction_quotient
    print(f'实际反应商 Q = {Q_actual:.2e}')
    print(f'实际驱动力 ΔG = {delta_g_actual:.2f} kJ/mol')
    
    # 判断反应方向
    if delta_g_actual < 0:
        print('在这些条件下，反应倾向于从左到右进行')
    elif delta_g_actual > 0:
        print('在这些条件下，反应倾向于从右到左进行')
    else:
        print('在这些条件下，反应处于平衡状态')

def example_with_custom_reaction():
    """
    使用自定义反应的示例
    """
    print('\\n' + '='*60)
    print('使用自定义代谢物浓度的示例')
    print('='*60)
    
    cc = ComponentContribution()
    
    # 一个典型的代谢反应: 葡萄糖 + ATP -> 葡萄糖-6-磷酸 + ADP
    rxn_str = "kegg:C00031 + kegg:C00002 = kegg:C00092 + kegg:C00008"  # glucose + ATP -> glucose-6-P + ADP
    reaction = cc.parse_reaction_formula(rxn_str)
    
    print(f'反应: {rxn_str}')
    
    # 生理条件下典型的代谢物浓度
    phys_conc = {
        'glucose': 5e-3,        # 5 mM
        'ATP': 1e-3,           # 1 mM
        'glucose-6-phosphate': 0.1e-3,  # 0.1 mM
        'ADP': 0.1e-3          # 0.1 mM
    }
    
    # 创建浓度字典（基于compound对象）
    concentrations = {}
    for compound, coeff in reaction.sparse.items():
        # 获取化合物名称
        compound_name = "Unknown"
        for identifier in compound.identifiers:
            if identifier.registry.namespace in ['kegg', 'bigg.metabolite', 'chebi', 'hmdb']:
                compound_name = identifier.accession.lower()
                break

        name_lower = compound_name.lower()
        if 'glucose' in name_lower and '6' not in name_lower and 'g6p' not in name_lower:
            concentrations[compound] = phys_conc['glucose']
        elif 'atp' in name_lower:
            concentrations[compound] = phys_conc['ATP']
        elif 'g6p' in name_lower or ('glucose' in name_lower and '6' in name_lower):
            concentrations[compound] = phys_conc['glucose-6-phosphate']
        elif 'adp' in name_lower:
            concentrations[compound] = phys_conc['ADP']
        else:
            concentrations[compound] = 1e-3  # 默认浓度
    
    print('\\n生理条件下的浓度:')
    for compound, conc in concentrations.items():
        # 获取化合物名称
        compound_name = "Unknown"
        for identifier in compound.identifiers:
            if identifier.registry.namespace in ['kegg', 'bigg.metabolite', 'chebi', 'hmdb']:
                compound_name = identifier.accession
                break
        print(f'  {compound_name}: {conc*1000:.3f} mM')
    
    # 设置条件并计算
    cc.p_h = Q_(7.0)
    cc.p_mg = Q_(3.0)
    cc.ionic_strength = Q_('0.25M')
    cc.temperature = Q_('298.15K')
    
    dg_standard = cc.standard_dg_prime(reaction).value.m_as('kJ/mol')
    
    # 计算浓度校正
    lnQ = 0
    for compound, coeff in reaction.sparse.items():
        lnQ += coeff * np.log(concentrations[compound])
    
    RT = R * 1e-3 * 298.15
    dg_actual = dg_standard + RT * lnQ
    
    print(f'\\n标准ΔG\'°: {dg_standard:.2f} kJ/mol')
    print(f'实际ΔG (生理浓度): {dg_actual:.2f} kJ/mol')
    print(f'浓度校正值: {RT * lnQ:.2f} kJ/mol')

if __name__ == '__main__':
    example_reaction_with_different_concentrations()
    example_with_custom_reaction()