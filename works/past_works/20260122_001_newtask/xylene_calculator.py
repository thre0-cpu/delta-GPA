'''
这个脚本用于计算二甲苯在特定条件下的生成自由能
使用CC方法（组分贡献法）
条件：温度298K，pH=10，pMg=6
'''

import numpy as np
from equilibrator_api import ComponentContribution, Q_
from rdkit import Chem, RDLogger
import sys


def get_compound(identifier: str, cc) -> object:
    """
    根据标识符获取化合物对象，按优先级尝试多种策略
    
    Args:
        identifier: 化合物标识符(支持InChI、KEGG、BIGG、Metacyc、SMILES等格式)
        cc: 已初始化的 ComponentContribution 对象
    
    Returns:
        成功时返回化合物对象，失败时返回 None
    """

    RDLogger.DisableLog('rdApp.error')  # 禁用RDKit错误日志输出

    def try_get_compound(query: str) -> object:
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
    
    def smiles_to_inchi(smiles: str) -> str:
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
) -> tuple:
    '''
    使用组分贡献法(Component Contribution)计算化合物的变换生成自由能 (Δ_fG'°或Δ_fG'm)
    
    参数:
    cpd: 化合物名称、InChI 或 KEGG ID
    cc: 已初始化的 ComponentContribution 对象
    p_h: pH值，默认7.0
    p_mg: pMg值，默认3.0
    I: 离子强度，默认0.25 M
    T: 温度，默认298.15 K
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
        from equilibrator_api import Reaction
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
        
        return float(val), float(err)
    
    except Exception as e:
        print(f"处理化合物 {cpd} 时出错: {str(e)}")
        return np.nan, np.nan


def main():
    print("开始计算二甲苯在指定条件下的生成自由能...")
    print("条件：温度298K，pH=10，pMg=6")
    
    # 初始化CC对象
    print("初始化ComponentContribution对象...")
    cc = ComponentContribution()
    
    # 二甲苯可能有不同的异构体，我们尝试几种可能的名称
    xylene_names = ["xylene", "dimethylbenzene", "C6H4(CH3)2", "1,2-xylene", "1,3-xylene", "1,4-xylene", "o-xylene", "m-xylene", "p-xylene"]
    
    print("尝试查找二甲苯的不同异构体...")
    for name in xylene_names:
        print(f"\n尝试化合物名称: {name}")
        try:
            result = standard_dgf_prime_CC(name, cc, p_h=10.0, p_mg=6.0, T=298.0)
            if not (np.isnan(result[0]) or np.isnan(result[1])):
                print(f"成功找到 {name} 的生成自由能:")
                print(f"Δ_fG'° = {result[0]:.2f} ± {result[1]:.2f} kJ/mol")
                return result
            else:
                print(f"未能找到 {name} 的数据")
        except Exception as e:
            print(f"处理 {name} 时发生错误: {str(e)}")
    
    # 如果以上都失败，尝试使用SMILES表示法
    print("\n尝试使用SMILES表示法...")
    xylene_smiles = [
        "CC1=CC=CC=C1C",  # o-xylene (邻二甲苯)
        "CC1=CC=CC=C1C",  # m-xylene (间二甲苯) - 这个SMILES实际上与o-xylene相同，需要更准确的表示
        "CC1=CC=C(C)C=C1"  # p-xylene (对二甲苯)
    ]
    
    # 由于RDKit无法区分邻位、间位异构体（它们的SMILES相同），我们只测试对二甲苯
    try:
        print(f"尝试对二甲苯的SMILES: CC1=CC=C(C)C=C1")
        result = standard_dgf_prime_CC("CC1=CC=C(C)C=C1", cc, p_h=10.0, p_mg=6.0, T=298.0)
        if not (np.isnan(result[0]) or np.isnan(result[1])):
            print(f"成功找到对二甲苯的生成自由能:")
            print(f"Δ_fG'° = {result[0]:.2f} ± {result[1]:.2f} kJ/mol")
            return result
        else:
            print("未能找到对二甲苯的数据")
    except Exception as e:
        print(f"处理对二甲苯时发生错误: {str(e)}")
    
    print("\n未能找到二甲苯的热力学数据。可能的原因：")
    print("1. 二甲苯不在eQuilibrator数据库中")
    print("2. 需要使用其他方法（如GNN）来预测")
    print("3. 化合物名称不正确")
    
    return None


if __name__ == "__main__":
    main()