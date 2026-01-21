import requests
import json
import numpy as np
from rdkit import Chem
from equilibrator_api import ComponentContribution, Q_
from typing import Optional, Tuple, Union
import pandas as pd
from tqdm import tqdm
import time

# 之前找到的含苯环化合物CID列表
benzene_cids = [980, 5143, 6643, 6758, 135398658, 5430, 243, 403, 785, 931, 991, 992, 996, 1486, 2345, 2912, 3314, 4685, 5794, 6294]

# 定义获取化合物属性的函数
def get_compound_properties(cids):
    """
    根据CID列表获取化合物的基本属性
    """
    if not cids:
        return []
    
    # 将CID列表转换为逗号分隔的字符串
    cid_str = ','.join(map(str, cids))
    
    # 构建属性请求URL
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_str}/property/MolecularFormula,InChI,InChIKey,CanonicalSMILES,IsomericSMILES,IUPACName,Title,XLogP,MolecularWeight/JSON"
    
    try:
        print("正在获取化合物属性...")
        response = requests.get(url)
        response.raise_for_status()
        
        data = response.json()
        properties = data.get('PropertyTable', {}).get('Properties', [])
        
        print(f"获取了 {len(properties)} 个化合物的属性")
        return properties
        
    except requests.exceptions.RequestException as e:
        print(f"获取化合物属性失败: {e}")
        return []
    except json.JSONDecodeError:
        print("响应不是有效的JSON格式")
        return []

def get_compound(identifier: str, cc):
    """
    根据标识符获取化合物对象，按优先级尝试多种策略
    
    Args:
        identifier: 化合物标识符
        cc: ChemicalCompound 或类似的化合物查询对象
    
    Returns:
        成功时返回化合物对象，失败时返回 None
    """
    def try_get_compound(query: str) -> Optional[object]:
        """尝试获取化合物，失败或返回None时返回None"""
        try:
            result = cc.get_compound(query)
            return result if result is not None else None
        except Exception:
            return None
    
    compound = None
    
    # 策略1: InChI
    if identifier.startswith("InChI="):
        try:
            compound = cc.get_compound_by_inchi(identifier)
        except Exception:
            return None
    
    # 策略2: KEGG
    elif identifier.startswith("C") and len(identifier) == 6 and identifier[1:].isdigit():
        compound = try_get_compound(f"kegg:{identifier}")
    
    # 策略3 & 4: BIGG 和 Metacyc
    if compound is None:
        # 尝试 BIGG
        compound = try_get_compound(f"bigg.metabolite:{identifier}")
        
        # BIGG 失败，尝试 Metacyc
        if compound is None:
            compound = try_get_compound(f"metacyc.compound:{identifier}")

            # Metacyc 失败，尝试搜索
            if compound is None:
                try:
                    compound = cc.search_compound(identifier)
                except Exception:
                    return None
    
    return compound

def standard_dgf_prime_CC(
    input: str, 
    p_h: float = 7.0, 
    p_mg: float = 3.0, 
    I: float = 0.25, 
    T: float = 298.15
) -> Tuple[np.floating, np.floating]:
    '''
    使用组分贡献法(Component Contribution)计算化合物的标准生成自由能

    注意:调用本函数需要同时调用 get_compound() 函数

    参数:
    input: 化合物的InChI字符串或其他Equilibrator API支持的格式
    p_h: 溶液的pH值 (默认值: 7.0)
    p_mg: 溶液的pMg值 (默认值: 3.0)
    I: 离子强度,单位为M (默认值: 0.25M)
    T: 温度,单位为K (默认值: 298.15K)
    
    返回:
    Tuple[np.floating, np.floating]: 
        - standard_dgf_prime_CC: 物质在指定条件下的生成自由能 (Δ_fG'°, kJ/mol)
        - std_CC: 生成自由能误差 (kJ/mol)
    '''
    cc = ComponentContribution()
    cc.p_h = Q_(p_h)
    cc.p_mg = Q_(p_mg)
    cc.ionic_strength = Q_(f'{I}M')
    cc.temperature = Q_(f'{T}K')

    # 获取化合物
    cpd = get_compound(input, cc)
    if cpd is None:
        raise ValueError(f"无法找到化合物: {input}")

    # 获取用户指定条件下的生化标准形成自由能 (Δ_fG'°)
    standard_dgf_prime_CC, sigma_fin, sigma_inf = cc.standard_dg_formation(cpd)
    
    # 使用 sigma_fin 作为有限不确定性
    std_CC = np.linalg.norm(sigma_fin) if sigma_fin is not None else np.float64(0.0)
    
    return np.float64(standard_dgf_prime_CC), np.float64(std_CC)

# 获取化合物详细属性
compound_properties = get_compound_properties(benzene_cids)

# 用CC方法计算pH=9下的生成自由能
def calculate_dgf_benzene_compounds(compound_properties, ph=9):
    """
    计算含有苯环的化合物的生成自由能
    
    Args:
        compound_properties: 从PubChem获取的化合物属性列表
        ph: pH值，默认为9
    
    Returns:
        包含化合物信息和计算结果的列表
    """
    results = []
    
    print(f"开始计算 {len(compound_properties)} 个化合物在pH={ph}下的生成自由能")
    
    for i, prop in enumerate(tqdm(compound_properties, desc="计算生成自由能")):
        cid = prop.get('CID')
        inchi = prop.get('InChI')
        name = prop.get('Title', 'Unknown')
        
        try:
            # 使用CC方法计算生成自由能
            # 注意：我们将pH设为9，其他参数使用默认值
            dgf, std = standard_dgf_prime_CC(inchi, p_h=ph, p_mg=3.0, I=0.25, T=298.15)
            
            result = {
                'CID': cid,
                'Name': name,
                'InChI': inchi,
                'SMILES': prop.get('CanonicalSMILES'),
                'MolecularWeight': prop.get('MolecularWeight'),
                'Dgf_prime': float(dgf),
                'Uncertainty': float(std),
                'pH': ph
            }
            
            results.append(result)
            
        except Exception as e:
            print(f"\n计算CID {cid} ({name}) 的生成自由能时出错: {e}")
            # 添加错误条目，便于后续分析
            result = {
                'CID': cid,
                'Name': name,
                'InChI': inchi,
                'SMILES': prop.get('CanonicalSMILES'),
                'MolecularWeight': prop.get('MolecularWeight'),
                'Dgf_prime': None,
                'Uncertainty': None,
                'pH': ph,
                'Error': str(e)
            }
            results.append(result)
    
    return results

# 执行计算
if compound_properties:
    print(f"\n获取到 {len(compound_properties)} 个化合物的属性")
    print("\n开始计算pH=9下的生成自由能")
    
    # 显示前几个化合物的信息
    print("\n前5个化合物的信息:")
    for i, prop in enumerate(compound_properties[:5]):
        print(f"{i+1}. CID: {prop.get('CID')} | 名称: {prop.get('Title', 'N/A')} | SMILES: {prop.get('CanonicalSMILES', 'N/A')}")
    
    calculation_results = calculate_dgf_benzene_compounds(compound_properties, ph=9)
    
    # 显示计算结果
    successful_calculations = [r for r in calculation_results if r['Dgf_prime'] is not None]
    print(f"\n成功计算了 {len(successful_calculations)} 个化合物的生成自由能")
    
    print("\npH=9下生成自由能计算结果 (成功计算的前10个):")
    for i, result in enumerate(successful_calculations[:10]):
        print(f"{i+1}. {result['Name']} (CID: {result['CID']})")
        print(f"   Dgf'°: {result['Dgf_prime']:.2f} kJ/mol ± {result['Uncertainty']:.2f}")
        print(f"   SMILES: {result['SMILES']}")
        print()
    
    # 保存结果到CSV文件
    if calculation_results:
        df = pd.DataFrame(calculation_results)
        df.to_csv('benzene_ring_compounds_dgf_pH9.csv', index=False)
        print("计算结果已保存到 benzene_ring_compounds_dgf_pH9.csv")
        
        # 显示结果摘要
        successful_df = df[df['Dgf_prime'].notna()]
        if not successful_df.empty:
            print("\n成功计算结果的摘要统计:")
            print(f"平均Dgf'°: {successful_df['Dgf_prime'].mean():.2f} kJ/mol")
            print(f"标准差: {successful_df['Dgf_prime'].std():.2f} kJ/mol")
            print(f"最小值: {successful_df['Dgf_prime'].min():.2f} kJ/mol")
            print(f"最大值: {successful_df['Dgf_prime'].max():.2f} kJ/mol")
        
        print(f"\n完整结果DataFrame的形状: {df.shape}")
        print("\n前几行结果:")
        print(df.head())