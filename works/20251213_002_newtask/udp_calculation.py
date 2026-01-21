import sys
import os
# 添加当前目录到Python路径
sys.path.append(os.getcwd())

# 导入必要的库
import numpy as np
from equilibrator_api import ComponentContribution, Q_
from equilibrator_cache import Reaction
from rdkit import Chem
from typing import Optional, Tuple, Union
import warnings

# 抑制可能的警告信息
warnings.filterwarnings('ignore')

def get_compound(identifier: str, cc) -> Optional[object]:
    """根据标识符获取化合物对象，按优先级尝试多种策略"""
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
    
    compound = None
    
    # 策略1: InChI
    if identifier.startswith('InChI='):
        try:
            compound = cc.get_compound_by_inchi(identifier)
        except Exception:
            return None
    
    # 策略2: KEGG
    elif identifier.startswith('C') and len(identifier) == 6 and identifier[1:].isdigit():
        compound = try_get_compound(f'kegg:{identifier}')
    
    # 策略3: SMILES
    elif is_smiles(identifier):
        try:
            inchi = smiles_to_inchi(identifier)
            if inchi:
                compound = cc.get_compound_by_inchi(inchi)
            else:
                return None
        except Exception:
            return None
    
    # 策略4 & 5: BIGG 和 Metacyc
    if compound is None:
        # 尝试 BIGG
        compound = try_get_compound(f'bigg.metabolite:{identifier}')
        
        # BIGG 失败，尝试 Metacyc
        if compound is None:
            compound = try_get_compound(f'metacyc.compound:{identifier}')

            # Metacyc 失败，尝试搜索
            if compound is None:
                try:
                    compound = cc.search_compound(identifier)
                except Exception:
                    return None
    
    return compound

def standard_dgf_prime_CC(
    input_str: str, 
    cc: ComponentContribution,
    p_h: float = 7.0, 
    p_mg: float = 3.0, 
    I: float = 0.25, 
    T: float = 298.15,
    physiological: bool = False # 是否转换为1mM生理浓度
) -> Tuple[np.floating, np.floating]:
    """
    使用组分贡献法(Component Contribution)计算化合物的变换生成自由能 (Δ_fG'°)
    
    参数:
    input_str: 化合物名称、InChI 或 KEGG ID
    cc: 已初始化的 ComponentContribution 对象
    p_h: pH 值 (默认值: 7.0)
    p_mg: pMg 值 (默认值: 3.0)
    I: 离子强度，单位为M (默认值: 0.25M)
    T: 温度，单位为K (默认值: 298.15K)
    physiological: 如果为 True，返回 1mM 浓度下的结果 (对应网页默认值)；
                   如果为 False (默认)，返回 1M 标准态结果。
    
    返回:
    (能量值, 误差值) 单位: kJ/mol
    """
    
    # 设置条件
    cc.p_h = Q_(p_h)
    cc.p_mg = Q_(p_mg)
    cc.ionic_strength = Q_(f'{I}M')
    cc.temperature = Q_(f'{T}K')

    # 获取化合物
    cpd = get_compound(input_str, cc)
    if cpd is None:
        raise ValueError(f'无法找到化合物: {input_str}')

    # 创建虚拟反应: 0 -> 1 化合物
    rxn = Reaction({cpd: 1})
    
    # 计算能量
    dg_prime_measurement = cc.standard_dg_prime(rxn)
    
    # 提取数值
    val = dg_prime_measurement.value.m  # 使用 .m 获取纯数值，避免单位问题
    err = dg_prime_measurement.error.m  # 使用 .m 获取纯数值，避免单位问题
    
    # 如果需要生理浓度 (1mM)，加上修正项 RT * ln(1e-3)
    if physiological:
        R = 8.314462618e-3  # kJ/(K·mol)
        correction = R * T * np.log(1e-3)
        val += correction
        # 误差通常不变，因为浓度修正项是确定值
    
    return np.float64(val), np.float64(err)

# 初始化 ComponentContribution 对象
# 注意: 这个过程可能比较慢，请耐心等待
print('正在初始化 ComponentContribution 对象...')
cc = ComponentContribution()
print('初始化完成！')

# 计算 UDP 在 298K, pH=10, pMg=6 条件下的生成自由能
print('\n开始计算 UDP 在 298K, pH=10, pMg=6 条件下的生成自由能...\n')

# 设置条件
temperature = 298.0  # K
ph = 10.0
pmg = 6.0
ionic_strength = 0.25  # M

# 计算UDP的生成自由能
try:
    dgf, err = standard_dgf_prime_CC(
        input_str='UDP',
        cc=cc,
        p_h=ph,
        p_mg=pmg,
        I=ionic_strength,
        T=temperature
    )
    
    print(f'UDP 在 {temperature}K, pH={ph}, pMg={pmg} 条件下的生成自由能:')
    print(f'Δ_fG\'° = {dgf:.2f} ± {err:.2f} kJ/mol')
    
    # 检查 UDP 化合物信息
    udp_compound = get_compound('UDP', cc)
    if udp_compound:
        print(f'\nUDP 化合物信息:')
        print(f'  ID: {udp_compound.id}')
        print(f'  化学式: {udp_compound.formula}')
        print(f'  名称: {udp_compound.name}')
        
except ValueError as e:
    print(f'错误: {e}')
    print('\n尝试使用其他标识符搜索UDP...')
    
    # 尝试其他可能的UDP标识符
    udp_identifiers = ['KEGG:C00015', 'bigg.metabolite:udp', 'Uridine diphosphate']
    
    for identifier in udp_identifiers:
        try:
            print(f'\n尝试使用标识符 \'{identifier}\'...')
            dgf, err = standard_dgf_prime_CC(
                input_str=identifier,
                cc=cc,
                p_h=ph,
                p_mg=pmg,
                I=ionic_strength,
                T=temperature
            )
            
            print(f'成功! UDP 在 {temperature}K, pH={ph}, pMg={pmg} 条件下的生成自由能:')
            print(f'Δ_fG\'° = {dgf:.2f} ± {err:.2f} kJ/mol')
            
            # 检查 UDP 化合物信息
            udp_compound = get_compound(identifier, cc)
            if udp_compound:
                print(f'\nUDP 化合物信息:')
                print(f'  ID: {udp_compound.id}')
                print(f'  化学式: {udp_compound.formula}')
                print(f'  名称: {udp_compound.name}')
            
            break
        
        except ValueError:
            print(f'标识符 \'{identifier}\' 未找到合适的UDP化合物')
            continue
    
    else:
        print('\n经过多种标识符尝试，仍未找到UDP化合物')