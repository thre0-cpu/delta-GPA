#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
使用CC方法计算生理浓度下的反应自由能变
"""

import pandas as pd
import numpy as np
from equilibrator_api import ComponentContribution, Q_
from rdkit import Chem, RDLogger
import warnings

warnings.filterwarnings('ignore')

def get_compound(identifier, cc):
    """
    根据标识符获取化合物对象，按优先级尝试多种策略
    
    Args:
        identifier: 化合物标识符(支持InChI、KEGG、BIGG、Metacyc、SMILES等格式)
        cc: 已初始化的 ComponentContribution 对象
    
    Returns:
        成功时返回化合物对象，失败时返回 None
    """
    RDLogger.DisableLog('rdApp.error')  # 禁用RDKit错误日志输出

    def try_get_compound(query):
        """尝试获取化合物，失败或返回None时返回None"""
        try:
            result = cc.get_compound(query)
            return result if result is not None else None
        except Exception:
            return None
    
    def is_smiles(s):
        """判断字符串是否为有效的SMILES格式"""
        try:
            mol = Chem.MolFromSmiles(s)
            return mol is not None
        except Exception:
            return False
    
    def smiles_to_inchi(smiles):
        """将SMILES转换为InChI"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                inchi = Chem.MolToInchi(mol)
                return inchi
            return None
        except Exception:
            return None
    
    # 策略0: 直接查询（原始标识符）
    compound = try_get_compound(identifier)
    if compound is not None:
        return compound
    
    # 策略0.5: 通用搜索
    try:
        compound = cc.search_compound(identifier)
        if compound is not None:
            return compound
    except Exception:
        pass

    # 策略1: InChI格式
    if identifier.startswith("InChI="):
        try:
            compound = cc.get_compound_by_inchi(identifier)
            if compound is not None:
                return compound
        except Exception:
            pass

    # 策略2: KEGG格式 (C + 5位数字)
    if identifier.startswith("C") and len(identifier) == 6 and identifier[1:].isdigit():
        compound = try_get_compound(f"kegg:{identifier}")
        if compound is not None:
            return compound
    
    # 策略3: BIGG数据库
    compound = try_get_compound(f"bigg.metabolite:{identifier}")
    if compound is not None:
        return compound
    
    # 策略4: MetaCyc数据库
    compound = try_get_compound(f"metacyc.compound:{identifier}")
    if compound is not None:
        return compound
    
    # 策略5: SMILES格式（最后尝试，因为检测成本较高）
    if is_smiles(identifier):
        try:
            inchi = smiles_to_inchi(identifier)
            if inchi:
                compound = try_get_compound(inchi)
                if compound is not None:
                    return compound
        except Exception:
            pass
    
    # 所有策略均失败
    return None

def standard_dgr_prime_CC_physiological(
    rxn: str,
    cc: ComponentContribution,
    p_h: float = 7.0,
    p_mg: float = 3.0,
    I: float = 0.25,
    T: float = 298.15,
    metabolite_conc: float = 1e-3  # 代谢物浓度，单位M，典型生理浓度是1mM
) -> tuple:
    """
    使用组分贡献法计算生理浓度下的反应自由能变

    参数:
    rxn: 反应式字符串
    cc: 已初始化的 ComponentContribution 对象
    p_h, p_mg, I, T: pH, pMg, 离子强度, 温度
    metabolite_conc: 代谢物浓度，单位M，默认为1mM (1e-3 M)

    返回:
    (能量值, 平衡常数) 单位: kJ/mol
    """
    try:
        # 设置条件
        cc.p_h = Q_(p_h)
        cc.p_mg = Q_(p_mg)
        cc.ionic_strength = Q_(f'{I}M')
        cc.temperature = Q_(f'{T}K')

        # 解析反应
        reaction = cc.parse_reaction_formula(rxn)

        # 计算标准反应的 ΔG'°
        dg_prime_measurement = cc.standard_dg_prime(reaction)

        # 提取数值
        dg_prime_std = dg_prime_measurement.value.m_as("kJ/mol")

        # 计算浓度校正项
        # ΔG = ΔG'° + RT ln(Q)，其中Q是反应商
        # 对于生理浓度 (1mM)，Q = (1e-3)^Σ(νi)，其中νi是反应系数
        R = 8.314462618e-3  # 气体常数，kJ/(mol*K)

        # 获取反应的化学计量系数
        # 在equilibrator_api中，使用reaction.sparse
        lnQ = 0
        for compound, coeff in reaction.sparse.items():
            # 浓度校正项是 coeff * ln([metabolite_conc])
            # coeff是化学计量系数，反应物为负，产物为正
            lnQ += coeff * np.log(metabolite_conc)

        # 应用浓度校正
        dg_prime_phys = dg_prime_std + R * T * lnQ

        # 计算平衡常数
        RT = R * T
        K_prime = np.exp(-dg_prime_std / RT)

        return float(dg_prime_phys), float(K_prime)

    except Exception as e:
        print(f"处理反应 {rxn} 时出错: {str(e)}")
        return np.nan, np.nan

def main():
    print('正在初始化ComponentContribution...')
    # 初始化CC对象
    cc = ComponentContribution()
    print('ComponentContribution初始化完成')
    
    # 读取数据文件
    df = pd.read_csv('TECRDB_test_filtered_full.csv')
    print(f'数据集包含 {len(df)} 个反应')
    
    # 计算生理条件下的反应自由能变
    print('开始计算生理浓度下的反应自由能变...')
    
    # 添加新列来存储生理条件下的计算结果
    df['physiological_dg_prime'] = np.nan
    df['physiological_k_prime'] = np.nan
    
    # 处理前几个反应进行演示
    for i in range(min(5, len(df))):
        print(f'处理第 {i+1} 个反应...')
        row = df.iloc[i]
        
        # 使用生理条件计算
        dg_phys, k_phys = standard_dgr_prime_CC_physiological(
            rxn=row['reaction'],
            cc=cc,
            p_h=row['p_h'],
            p_mg=row['p_mg'],
            I=row['ionic_strength'],
            T=row['temperature'],
            metabolite_conc=1e-3  # 1mM生理浓度
        )
        
        df.loc[i, 'physiological_dg_prime'] = dg_phys
        df.loc[i, 'physiological_k_prime'] = k_phys
        
        print(f'  反应: {row["reaction"]}')
        print(f'  标准ΔG\'°: {np.log(row["K_prime"]) * -8.314e-3 * row["temperature"]:.2f} kJ/mol')
        print(f'  标准K\': {row["K_prime"]:.2e}')
        print(f'  计算ΔG\'°: {np.log(row["calculated_K_prime"]) * -8.314e-3 * row["temperature"]:.2f} kJ/mol')
        print(f'  计算K\': {row["calculated_K_prime"]:.2e}')
        print(f'  生理ΔG\'°: {dg_phys:.2f} kJ/mol')
        print(f'  生理K\': {k_phys:.2e}')
        print()
    
    # 保存结果
    df.to_csv('TECRDB_test_with_physiological_predictions.csv', index=False)
    print('生理条件下的预测结果已保存到 TECRDB_test_with_physiological_predictions.csv')
    
    # 显示几个主要统计量的比较
    print('生理条件下计算结果统计:')
    print(f'生理ΔG\'°范围: {df["physiological_dg_prime"].min():.2f} - {df["physiological_dg_prime"].max():.2f} kJ/mol')
    print(f'生理K\'范围: {df["physiological_k_prime"].min():.2e} - {df["physiological_k_prime"].max():.2e}')

if __name__ == '__main__':
    main()