"""
苯环化合物生成自由能计算脚本

本脚本用于从PubChem数据库筛选含有苯环的化合物，并使用CC方法计算它们在pH=9下的生成自由能。
"""

import pubchempy as pcp
import pandas as pd
import numpy as np
import requests
from tqdm import tqdm
import time
import re
import json
import os
from typing import Optional

import numpy as np
from rdkit import Chem
from equilibrator_api import ComponentContribution, Q_

def get_compound(identifier: str, cc) -> Optional[object]:
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


def standard_dgf_prime_CC(input, p_h=7.0, p_mg=3.0, I=0.25, T=298.15):
    '''
    使用组分贡献法（Component Contribution）计算化合物的标准生成自由能

    注意：调用本函数需要同时调用 get_compound() 函数

    参数:
    input: 化合物的InChI字符串或其他Equilibrator API支持的格式
    p_h: 溶液的pH值 (默认值: 7.0)
    p_mg: 溶液的pMg值 (默认值: 3.0)
    I: 离子强度，单位为M (默认值: 0.25M)
    T: 温度，单位为K (默认值: 298.15K)
    
    返回:
    standard_dgf_prime_CC: 物质在指定条件下的生成自由能 (Δ_fG'°, kJ/mol)
    std_CC: 生成自由能误差 (kJ/mol)
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
    std_CC = np.linalg.norm(sigma_fin) if sigma_fin is not None else 0.0
    
    return standard_dgf_prime_CC, std_CC


def has_benzene_ring(smiles):
    """
    检查SMILES字符串是否含有苯环
    """
    # 检查是否含有苯环的模式
    # 这包括各种可能的表达方式，如c1ccccc1, c1ccc(cc1), C1=CC=CC=C1等
    benzene_patterns = [
        r'c1ccccc1',  # 经典的苯环
        r'c1ccc([cC]c1)',  # 苯环上有取代基
        r'[cC]1=[cC][cC]=[cC][cC]=[cC]1',  # 双键表示的苯环
    ]
    
    for pattern in benzene_patterns:
        if re.search(pattern, smiles, re.IGNORECASE):
            return True
    return False


def get_benzene_compounds(limit=100):
    """
    通过结构搜索获取含苯环的化合物
    """
    print(f"正在搜索含有苯环的化合物，限制数量为 {limit} ...")

    # 使用苯环作为子结构进行搜索
    # 我们使用SMILES:c1ccccc1来表示苯环
    # 但由于苯环是非常常见的结构，子结构搜索会返回大量结果，导致JSON解析错误
    # 我们需要限制返回的CIDs数量
    benzene_smiles = 'c1ccccc1'

    try:
        # 使用POST请求，限制返回的CIDs数量
        # 由于PubChem API的限制，直接限制CIDs数量可能不可行
        # 我们先获取一个小样本，再从样本中获取化合物信息
        headers = {
            'Content-Type': 'application/x-www-form-urlencoded'
        }

        # 尝试使用URL参数限制结果
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/smiles/{benzene_smiles}/property/Title,SMILES,InChI/JSON"

        # 由于直接限制数量可能困难，我们先获取一个大的结果集然后限制处理的数量
        response = requests.get(url)

        if response.status_code == 200:
            # 尝试解析响应，如果失败则说明数据量太大
            try:
                data = response.json()
                properties = data.get('PropertyTable', {}).get('Properties', [])
                valid_cids = []

                for prop in properties:
                    cid = prop.get('CID')
                    smiles = prop.get('SMILES', '')

                    # 验证该化合物确实包含苯环
                    if has_benzene_ring(smiles):
                        valid_cids.append(cid)
                        if len(valid_cids) >= limit:
                            break

                print(f"找到 {len(valid_cids)} 个含有苯环的化合物")
                return valid_cids

            except json.JSONDecodeError as e:
                print(f"JSON解析错误，可能是因为数据量过大: {e}")
                # 如果JSON解析失败，尝试获取CIDs列表，而不是完整的属性
                cids_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/smiles/{benzene_smiles}/cids/JSON"
                cids_response = requests.get(cids_url)

                if cids_response.status_code == 200:
                    cids_data = cids_response.json()
                    all_cids = cids_data.get('IdentifierList', {}).get('CID', [])

                    print(f"获取到 {len(all_cids)} 个CIDs，将限制到前 {limit} 个进行处理")

                    # 限制处理的数量以避免API限制
                    limited_cids = all_cids[:min(limit*5, len(all_cids))]
                    return limited_cids
                else:
                    print(f"获取CIDs列表失败，状态码: {cids_response.status_code}")
                    return []
        else:
            print(f"API请求失败，状态码: {response.status_code}")
            return []

    except Exception as e:
        print(f"获取苯环化合物时发生错误: {e}")
        return []


def get_compound_details(cids):
    """
    获取化合物的详细信息
    """
    print(f"正在获取 {len(cids)} 个化合物的详细信息...")
    details = []
    
    # 批量处理，避免请求过于频繁
    for i in tqdm(range(0, len(cids), 50), desc="获取化合物信息"):
        batch = cids[i:i+50]
        batch_cids_str = ','.join(map(str, batch))
        
        try:
            # 获取化合物的SMILES和名称
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{batch_cids_str}/property/Title,SMILES,InChI/JSON"
            response = requests.get(url)
            
            if response.status_code == 200:
                data = response.json()
                properties = data.get('PropertyTable', {}).get('Properties', [])
                
                for prop in properties:
                    cid = prop.get('CID')
                    title = prop.get('Title', 'N/A')
                    smiles = prop.get('SMILES', 'N/A')
                    inchi = prop.get('InChI', 'N/A')
                    
                    details.append({
                        'CID': cid,
                        'Title': title,
                        'SMILES': smiles,
                        'InChI': inchi,
                        'HasBenzeneRing': has_benzene_ring(smiles)
                    })
                    
            else:
                print(f"获取化合物信息失败，状态码: {response.status_code}, CIDs: {batch_cids_str}")
            
            # 限制请求频率，避免超过API限制
            time.sleep(0.5)
            
        except Exception as e:
            print(f"获取化合物信息时发生错误: {e}, CIDs: {batch_cids_str}")
            continue
    
    return details


def calculate_dgf_with_CC(compound_info, pH=9):
    """
    使用CC方法计算化合物在指定pH下的生成自由能
    """
    try:
        # 获取化合物的InChI表示
        inchi = compound_info['InChI']
        
        # 调用CC方法计算生成自由能
        dgf_prime, uncertainty = standard_dgf_prime_CC(inchi, p_h=pH, p_mg=3.0, I=0.25, T=298.15)
        
        return dgf_prime.m, uncertainty
    except Exception as e:
        print(f"计算 {compound_info['Title']} ({compound_info['CID']}) 的生成自由能时出错: {e}")
        return None, None


def main():
    # 获取含苯环的化合物ID列表
    benzene_compound_cids = get_benzene_compounds(limit=100)
    print(f"获取了 {len(benzene_compound_cids)} 个苯环化合物用于后续计算")
    
    if not benzene_compound_cids:
        print("未能获取任何含苯环的化合物，程序退出")
        return

    # 获取化合物详细信息
    compound_details = get_compound_details(benzene_compound_cids)

    # 过滤出确认含有苯环的化合物
    benzene_compounds = [comp for comp in compound_details if comp['HasBenzeneRing']]
    print(f"确认含有苯环的化合物数量: {len(benzene_compounds)}")

    # 创建DataFrame并保存
    df = pd.DataFrame(benzene_compounds)
    df.to_csv('benzene_compounds.csv', index=False)
    print(f"已保存 {len(df)} 个含苯环化合物的信息到 benzene_compounds.csv")

    # 计算所有苯环化合物的生成自由能
    results = []

    for i, compound in enumerate(tqdm(benzene_compounds, desc="计算生成自由能 (pH=9)")):
        dgf, uncertainty = calculate_dgf_with_CC(compound, pH=9)
        
        results.append({
            'CID': compound['CID'],
            'Title': compound['Title'],
            'SMILES': compound['SMILES'],
            'InChI': compound['InChI'],
            'DGF_prime_pH9': dgf,
            'Uncertainty': uncertainty
        })
        
        # 显示进度
        if (i+1) % 10 == 0:
            print(f"已处理 {i+1}/{len(benzene_compounds)} 个化合物")

    # 创建结果DataFrame
    results_df = pd.DataFrame(results)

    # 移除计算失败的条目
    results_df = results_df.dropna(subset=['DGF_prime_pH9'])
    print(f"成功计算了 {len(results_df)} 个化合物的生成自由能")

    # 保存结果
    results_df.to_csv('benzene_compounds_dgf_pH9.csv', index=False)
    print("结果已保存到 benzene_compounds_dgf_pH9.csv")

    # 显示前几行结果
    print("\n前10个化合物的生成自由能 (pH=9):")
    print(results_df[['CID', 'Title', 'SMILES', 'DGF_prime_pH9', 'Uncertainty']].head(10))

    # 基本统计信息
    print("\n生成自由能 (pH=9) 统计信息:")
    if 'DGF_prime_pH9' in results_df.columns:
        print(f"平均值: {results_df['DGF_prime_pH9'].mean():.2f} kJ/mol")
        print(f"标准差: {results_df['DGF_prime_pH9'].std():.2f} kJ/mol")
        print(f"最小值: {results_df['DGF_prime_pH9'].min():.2f} kJ/mol")
        print(f"最大值: {results_df['DGF_prime_pH9'].max():.2f} kJ/mol")
        print(f"数量: {len(results_df)}")
    else:
        print("未能计算任何化合物的生成自由能")


if __name__ == "__main__":
    main()