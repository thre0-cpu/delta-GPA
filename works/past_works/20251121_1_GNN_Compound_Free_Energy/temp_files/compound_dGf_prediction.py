'''
使用GNN模型预测化合物 'CCCCCCCC(O)CCCCCC(O)CCCC(=O)O' 的标准生成自由能
'''
import subprocess
import json
import sys
import os
from typing import Optional, Tuple, Union, List

import numpy as np
import numpy.typing as npt
import pandas as pd
from rdkit import Chem

from equilibrator_api import ComponentContribution, Q_

from rdkit import RDLogger
import warnings

# 抑制报错输出
RDLogger.DisableLog('rdApp.*')
warnings.filterwarnings('ignore', category=DeprecationWarning)

# 初始化CC类，以便在后续函数中调用
cc = ComponentContribution() # 初始化组分贡献法对象放在全局，避免重复初始化

# 导入GNN预测函数
from APIs.Prediction_tools.GNN.model.GNN_dGf import predict_standard_dGf_prime


def convert_to_inchi(input_string):
    """
    输入字符串，如果是SMILES就转换成InChI，如果是InChI就直接保留
    
    参数:
        input_string: 输入的化学结构字符串
    
    返回:
        InChI格式的字符串
    """
    input_string = input_string.strip()
    
    # 判断是否为InChI格式（InChI以"InChI="开头）
    if input_string.startswith("InChI="):
        print("检测到InChI格式，直接保留")
        return input_string
    else:
        # 假设是SMILES格式，尝试转换为InChI
        print("检测到SMILES格式，转换为InChI")
        try:
            mol = Chem.MolFromSmiles(input_string, sanitize=False)
            if mol is not None:
                inchi = Chem.MolToInchi(mol)
                return inchi
            else:
                return "错误：无效的SMILES字符串"
        except Exception as e:
            return f"错误：转换失败 - {str(e)}"


def get_compound(identifier: str, cc) -> Optional[object]:
    """
    根据标识符获取化合物对象，按优先级尝试多种策略
    
    优先级顺序：
    1. InChI 格式
    2. SMILES 格式（转换为 InChI）
    3. KEGG ID
    4. BIGG ID
    5. Metacyc ID
    6. 名称搜索
    """
    def try_get_compound(query: str) -> Optional[object]:
        try:
            result = cc.get_compound(query)
            return result if result is not None else None
        except Exception:
            return None
    
    def is_smiles(s: str) -> bool:
        try:
            mol = Chem.MolFromSmiles(s, sanitize=False)
            return mol is not None
        except Exception:
            return False
    
    def smiles_to_inchi(smiles: str) -> Optional[str]:
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is not None:
                return Chem.MolToInchi(mol)
            return None
        except Exception:
            return None
    
    # 策略1: InChI 格式
    if identifier.startswith("InChI="):
        try:
            compound = cc.get_compound_by_inchi(identifier)
            if compound:
                return compound
        except Exception:
            pass
    
    # 策略2: SMILES 格式（优先处理结构化数据）
    if is_smiles(identifier):
        try:
            inchi = smiles_to_inchi(identifier)
            if inchi:
                compound = cc.get_compound_by_inchi(inchi)
                if compound:
                    return compound
        except Exception:
            pass
    
    # 策略3: KEGG ID
    if identifier.startswith("C") and len(identifier) == 6 and identifier[1:].isdigit():
        compound = try_get_compound(f"kegg:{identifier}")
        if compound:
            return compound
    
    # 策略4: BIGG ID
    compound = try_get_compound(f"bigg.metabolite:{identifier}")
    if compound:
        return compound
    
    # 策略5: Metacyc ID
    compound = try_get_compound(f"metacyc.compound:{identifier}")
    if compound:
        return compound
    
    # 策略6: 名称搜索（最后尝试）
    try:
        compound = cc.search_compound(identifier)
        if compound:
            return compound
    except Exception:
        pass
    
    return None


def calculate_standard_dgf_GNN(
    identifier: str,
    cc,
    p_h: float = 7.0,
    p_mg: float = 3.0,
    I: float = 0.25,
    T: float = 298.15
) -> Tuple[np.floating, np.floating, str]:
    """
    自动判断化合物是否已知，并选择合适的方法计算标准生成自由能
    
    工作流程:
    1. 尝试通过 get_compound() 查询化合物
    2. 如果找到(已知化合物)：使用 cc_transformed_standard_dgf_prime_GNN (经过Legendre变换)
    3. 如果未找到(未知化合物)：使用 untransformed_standard_dgf_GNN (未经Legendre变换)
    
    注意：调用本函数必须先全局初始化cc = ComponentContribution()，同时调用以下函数：
        get_compound()
        cc_transformed_standard_dgf_prime_GNN()
        convert_to_inchi()
        untransformed_standard_dgf_GNN()

    Args:
        identifier: 化合物标识符（支持InChI、KEGG、BIGG、Metacyc、SMILES等格式）
        cc: ComponentContribution 对象
        p_h: pH值，默认7.0
        p_mg: pMg值，默认3.0
        I: 离子强度，默认0.25 M
        T: 温度，默认298.15 K
    
    Returns:
        Tuple[np.floating, np.floating, str]: 
            - 标准生成自由能 (kJ/mol)
            - 不确定度 (kJ/mol)
            - 计算方法标识 ("known_transformed" 或 "unknown_untransformed")
    
    Example:
        >>> cc = ComponentContribution()
        >>> # 已知化合物
        >>> dgf, uncertainty, method = calculate_standard_dgf_auto("KEGG:C00002", cc)
        >>> print(f"Method: {method}, ΔGf: {dgf:.2f} ± {uncertainty:.2f} kJ/mol")
        
        >>> # 未知化合物
        >>> smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"  # 布洛芬
        >>> dgf, uncertainty, method = calculate_standard_dgf_auto(smiles, cc)
        >>> print(f"Method: {method}, ΔGf: {dgf:.2f} ± {uncertainty:.2f} kJ/mol")
    """
    # 尝试获取化合物对象
    compound = get_compound(identifier, cc)
    
    if compound is not None:
        # 检查化合物是否有有效的 InChI
        has_valid_inchi = hasattr(compound, 'inchi') and compound.inchi is not None and compound.inchi.strip() != ""
        
        if has_valid_inchi:
            # 已知化合物且有InChI：使用经过Legendre变换的方法（CC数据库）
            try:
                dgf, uncertainty = cc_transformed_standard_dgf_prime_GNN(
                    input=identifier,
                    cc=cc,
                    p_h=p_h,
                    p_mg=p_mg,
                    I=I,
                    T=T
                )
                method = "known_transformed"
                print(f"✓ 已知化合物：使用CC数据库的Legendre变换方法")
                return dgf, uncertainty, method
                
            except Exception as e:
                print(f"⚠ 已知化合物但计算失败: {str(e)}")
                print(f"  尝试作为未知化合物处理...")
                # 继续执行下面的未知化合物逻辑
        else:
            print(f"⚠ 化合物已找到但缺少InChI信息")
            print(f"  尝试作为未知化合物处理...")
    
    # 未知化合物或已知化合物但缺少InChI：使用未经Legendre变换的GNN预测方法
    print(f"✗ 未知化合物：使用GNN预测方法（未经Legendre变换，缺少pKa信息）")
    
    # 需要先转换为InChI格式
    inchi = convert_to_inchi(identifier)
    
    dgf, uncertainty = untransformed_standard_dgf_GNN(
        input=inchi,
        cc=cc,
        p_h=p_h,
        p_mg=p_mg,
        I=I,
        T=T
    )
    method = "unknown_untransformed"
    
    return dgf, uncertainty, method


def cc_transformed_standard_dgf_prime_GNN(
    input: str, 
    cc,
    p_h: float = 7.0, 
    p_mg: float = 3.0, 
    I: float = 0.25, 
    T: float = 298.15
) -> Tuple[np.floating, np.floating]:
    '''
    使用图神经网络(GNN)计算已知化合物的标准生成自由能（经过CC数据库的Legendre变换）
    
    注意：调用本函数必须先全局初始化cc = ComponentContribution()，同时调用 get_compound() 函数
    
    参数:
    input: 化合物的InChI字符串或其他Equilibrator API支持的格式
    p_h: 溶液的pH值 (默认值: 7.0)
    p_mg: 溶液的pMg值 (默认值: 3.0)
    I: 离子强度，单位为M (默认值: 0.25M)
    T: 温度，单位为K (默认值: 298.15K)
    
    返回:
    standard_dgf_prime_GNN: 物质在指定条件下的形成自由能 (Δ_fG'°, kJ/mol)
    std_GNN: 生成自由能误差 (kJ/mol)
    '''
    
    # 获取化合物
    cpd = get_compound(input, cc)
    if cpd is None:
        raise ValueError(f"无法找到化合物: {input}")

    # 步骤1: GNN 返回 未 transform 的自由能值
    standard_dgf_c_GNN, std_GNN = predict_standard_dGf_prime(cpd.inchi)
    
    # 步骤2: Legendre变换，转换到用户指定的条件
    delta_user_c = cpd.transform(p_h = Q_(7.0), p_mg = Q_(14.0), ionic_strength = Q_('0.25M'), temperature = Q_('298.15K'))
    standard_dgf_GNN = standard_dgf_c_GNN - delta_user_c.m_as("kJ/mol")

    delta_user = cpd.transform(p_h = Q_(p_h), p_mg = Q_(p_mg), ionic_strength = Q_(f'{I}M'), temperature = Q_(f'{T}K'))
    standard_dgf_prime_GNN = standard_dgf_GNN + delta_user.m_as("kJ/mol")

    return float(standard_dgf_prime_GNN), float(std_GNN)


def untransformed_standard_dgf_GNN(
    input: str, 
    cc,
    p_h: float = 7.0, 
    p_mg: float = 3.0, 
    I: float = 0.25, 
    T: float = 298.15
) -> Tuple[np.floating, np.floating]:
    '''
    使用图神经网络(GNN)计算全新化合物的生成自由能（缺少pKa信息，未经过Legendre变换）
    
    注意：调用本函数必须先全局初始化cc = ComponentContribution()，同时调用 convert_to_inchi() 函数
    
    参数:
    input: 化合物的InChI字符串或其他Equilibrator API支持的格式
    p_h: 溶液的pH值 (默认值: 7.0)
    p_mg: 溶液的pMg值 (默认值: 3.0)
    I: 离子强度，单位为M (默认值: 0.25M)
    T: 温度，单位为K (默认值: 298.15K)
    
    返回:
    standard_dgf_prime_GNN: 物质在指定条件下的形成自由能 (Δ_fG'°, kJ/mol)
    std_GNN: 生成自由能误差 (kJ/mol)
    '''

    cpd_inchi = convert_to_inchi(input)
    # GNN 返回 未 transform 的自由能值
    standard_dgf_GNN, std_GNN = predict_standard_dGf_prime(cpd_inchi)

    return float(standard_dgf_GNN), float(std_GNN)


# 主程序
if __name__ == "__main__":
    print("开始计算化合物 'CCCCCCCC(O)CCCCCC(O)CCCC(=O)O' 的标准生成自由能")
    print("="*60)
    
    # 化合物标识符
    compound_identifier = 'CCCCCCCC(O)CCCCCC(O)CCCC(=O)O'
    
    # 首先尝试将SMILES转换为InChI
    print("1. 将SMILES转换为InChI:")
    inchi = convert_to_inchi(compound_identifier)
    print(f"   输入: {compound_identifier}")
    print(f"   InChI: {inchi}")
    print()
    
    # 使用GNN模型计算标准生成自由能
    print("2. 使用GNN模型计算标准生成自由能:")
    try:
        dgf, uncertainty, method = calculate_standard_dgf_GNN(compound_identifier, cc)
        print()
        print("="*60)
        print("结果:")
        print(f"   化合物: {compound_identifier}")
        print(f"   ΔGf'°: {dgf:.2f} kJ/mol")
        print(f"   不确定度: ±{uncertainty:.2f} kJ/mol")
        print(f"   计算方法: {method}")
        print("="*60)
    except Exception as e:
        print(f"计算过程中出现错误: {str(e)}")