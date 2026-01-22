#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
使用CC方法预测TECRDB_test.csv中反应的平衡常数
本脚本用于使用Component Contribution (CC)方法预测TECRDB_test.csv文件中反应的平衡常数。
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


def calculate_equilibrium_constant(row, cc):
    """
    使用CC方法计算单个反应的平衡常数
    """
    try:
        # 解析反应
        reaction = cc.parse_reaction_formula(row['reaction'])

        # 保存原始条件
        original_pH = cc.p_h
        original_pMg = cc.p_mg
        original_ionic_strength = cc.ionic_strength
        original_temperature = cc.temperature

        # 设置条件
        cc.p_h = Q_(row['p_h'])
        cc.p_mg = Q_(row['p_mg'])
        cc.ionic_strength = Q_(f"{row['ionic_strength']}M")
        cc.temperature = Q_(f"{row['temperature']}K")

        # 计算反应的 ΔG'°
        dg_prime = cc.standard_dg_prime(reaction)

        # 恢复原始条件
        cc.p_h = original_pH
        cc.p_mg = original_pMg
        cc.ionic_strength = original_ionic_strength
        cc.temperature = original_temperature

        # 获取数值（注意：返回的是Measurement对象，需要提取值）
        dg_prime_kJ_per_mol = dg_prime.value.m_as('kJ/mol')

        # 计算平衡常数
        # ΔG'° = -RT ln(K')
        # K' = exp(-ΔG'°/RT)
        RT = 8.314462618e-3 * row['temperature']  # R*T, R = 8.314462618 J/(mol*K)，转换为kJ/(mol*K)
        K_prime = np.exp(-dg_prime_kJ_per_mol / RT)

        return dg_prime_kJ_per_mol, K_prime

    except Exception as e:
        print(f"处理反应 {row['reaction']} 时出错: {str(e)}")
        return np.nan, np.nan


def main():
    print('正在初始化ComponentContribution...')
    # 初始化CC对象
    cc = ComponentContribution()
    print('ComponentContribution初始化完成')

    # 读取数据文件
    df = pd.read_csv('TECRDB_test.csv')
    print(f'数据集包含 {len(df)} 个反应')
    print('数据集的列包括：', list(df.columns))
    print('\n前5行数据：')
    print(df.head())

    print('开始计算每个反应的平衡常数...')
    
    # 添加新列来存储计算结果
    df['calculated_dg_prime'] = np.nan
    df['calculated_K_prime'] = np.nan

    # 处理前几行进行测试
    for i in range(min(10, len(df))):
        print(f'处理第 {i+1} 行反应...')
        row = df.iloc[i]
        dg_prime, K_prime = calculate_equilibrium_constant(row, cc)
        
        df.loc[i, 'calculated_dg_prime'] = dg_prime
        df.loc[i, 'calculated_K_prime'] = K_prime
        
        print(f'  反应: {row["reaction"]}')
        print(f'  计算的ΔG\'°: {dg_prime:.2f} kJ/mol')
        print(f'  计算的K\': {K_prime:.2e}')
        print(f'  文献K\': {row["K_prime"]}')
        print()

    print('前10个反应计算完成')

    # 保存结果
    df.to_csv('TECRDB_test_with_predictions.csv', index=False)
    print('结果已保存到 TECRDB_test_with_predictions.csv')


if __name__ == '__main__':
    main()