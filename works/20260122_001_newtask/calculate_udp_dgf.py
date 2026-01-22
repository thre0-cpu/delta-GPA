'''
计算UDP在298K，pH=10，pMg=6条件下的生成自由能
使用组分贡献法(Component Contribution)
'''

import subprocess
import json
import sys
import os
from typing import Optional, Tuple, Union

import numpy as np
import numpy.typing as npt
from rdkit import Chem, RDLogger

from equilibrator_api import ComponentContribution, Reaction, Q_

def get_compound(identifier: str, cc) -> Optional[object]:
    """
    根据标识符获取化合物对象，按优先级尝试多种策略

    注意：调用本函数必须全局初始化cc = ComponentContribution()
    
    Args:
        identifier: 化合物标识符(支持InChI、KEGG、BIGG、Metacyc、SMILES等格式)
        - 格式1: InChI (以 "InChI=" 开头)
        - 格式2: KEGG (以 "C" 开头，后跟5位数字)
        - 格式3: BIGG (直接使用BIGG ID)
        - 格式4: Metacyc (直接使用Metacyc ID)
        - 格式5: SMILES (有效的SMILES字符串，最后尝试)
        cc: 已初始化的 ComponentContribution 对象
    
    Returns:
        成功时返回化合物对象，失败时返回 None
    """

    RDLogger.DisableLog('rdApp.error')  # 禁用RDKit错误日志输出

    def try_get_compound(query: str) -> Optional[object]:
        """尝试获取化合物，失败或返回None时返回None"""
        try:
            result = cc.get_compound(query)
            return result if result is not None else None
        except Exception:
            return None
    
    def is_smiles(s: str) -> bool:
        """判断字符串是否为有效的SMILES格式"""
        try:
            mol = Chem.MolFromSmiles(s)
            return mol is not None
        except Exception:
            return False
    
    def smiles_to_inchi(smiles: str) -> Optional[str]:
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
    compound = cc.search_compound(identifier)
    if compound is not None:
        return compound

    
    # 策略1: InChI格式
    compound = cc.get_compound_by_inchi(identifier)
    if compound is not None:
        return compound


    # 策略2: KEGG格式 (C + 5位数字)
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
                compound = cc.get_compound_by_inchi(inchi)
                if compound is not None:
                    return compound
        except Exception:
            pass
    
    # 所有策略均失败
    return None

def standard_dgf_prime_CC(
    cpd: str, 
    cc: ComponentContribution,
    p_h: float = 7.0, 
    p_mg: float = 3.0, 
    I: float = 0.25, 
    T: float = 298.15,
    physiological: bool = False # 是否转换为1mM生理浓度
) -> Tuple[np.floating, np.floating]:
    '''
    使用组分贡献法(Component Contribution)计算化合物的变换生成自由能 (Δ_fG'°或Δ_fG'm)

    注意：调用本函数必须提前全局初始化cc = ComponentContribution()，并调用get_compound()获取化合物对象。
    
    参数:
    cpd: 化合物名称、InChI 或 KEGG ID
    cc: 已初始化的 ComponentContribution 对象
    physiological: 如果为 True，返回 1mM 浓度下的结果；
                   如果为 False (默认)，返回 1M 标准态结果。
    
    返回:
    (能量值, 误差值) 单位: kJ/mol
    '''
    try:
        # 设置条件
        cc.p_h = Q_(p_h)
        cc.p_mg = Q_(p_mg)
        cc.ionic_strength = Q_(f'{I}M')
        cc.temperature = Q_(f'{T}K')

        # 获取化合物
        compound = get_compound(cpd, cc)
        if compound is None:
            raise ValueError(f"无法找到化合物: {cpd}")

        # 创建虚拟反应: 0 -> 1 化合物
        rxn_c = Reaction({compound: 1})
        
        # 计算能量
        dg_prime_measurement = cc.standard_dg_prime(rxn_c)
        
        # 提取数值
        val = dg_prime_measurement.value.m_as("kJ/mol")
        err = dg_prime_measurement.error.m_as("kJ/mol")
        
        # 如果需要生理浓度 (1mM)，加上修正项
        if physiological:
            from scipy.constants import R # J/(K·mol)
            R = R * 1e-3  # kJ/(K·mol)
            correction = R * T * np.log(1e-3)
            val += correction
        
        return np.float64(val), np.float64(err)
    
    except Exception as e:
        print(f"处理化合物 {cpd} 时出错: {str(e)}")
        return np.nan, np.nan

# 初始化CC类
print("初始化ComponentContribution...")
cc = ComponentContribution()

# 计算UDP在298K，pH=10，pMg=6条件下的生成自由能
print("计算UDP在298K，pH=10，pMg=6条件下的生成自由能")
dGf, std = standard_dgf_prime_CC('UDP', cc, p_h=10, p_mg=6, T=298.15)
print(f"UDP Δ_fG'° (pH=10, pMg=6, T=298.15K): {dGf:.2f} ± {std:.2f} kJ/mol")