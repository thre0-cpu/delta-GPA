import subprocess
import json
import sys
import os
import requests
import time
import pandas as pd
import numpy as np
from typing import Optional
from tqdm import tqdm

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


def standard_dgf_prime_CC(input, cc, p_h=7.0, p_mg=3.0, I=0.25, T=298.15):
    """
    使用组分贡献法（Component Contribution）计算化合物的标准生成自由能

    注意：调用本函数需要同时调用 get_compound() 函数

    参数:
    input: 化合物的InChI字符串或其他Equilibrator API支持的格式
    cc: ComponentContribution 实例，用于避免重复初始化
    p_h: 溶液的pH值 (默认值: 7.0)
    p_mg: 溶液的pMg值 (默认值: 3.0)
    I: 离子强度，单位为M (默认值: 0.25M)
    T: 温度，单位为K (默认值: 298.15K)
    
    返回:
    standard_dgf_prime_CC: 物质在指定条件下的生成自由能 (Δ_fG'°, kJ/mol)
    std_CC: 生成自由能误差 (kJ/mol)
    """
    # 设置条件
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


# 使用PubChem API搜索含苯环的化合物
def search_benzene_compounds(limit=20):
    """搜索含苯环的化合物"""
    # 使用POST方法进行子结构搜索，因为SMILES包含特殊字符
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/substructure/smiles/cid/TXT"
    smiles = "c1ccccc1"  # 苯环SMILES
    payload = smiles

    # 设置Content-Type header
    headers = {
        'Content-Type': 'text/plain'
    }

    response = requests.post(url, data=payload, headers=headers)

    if response.status_code == 200:
        # 响应是纯文本格式的CID列表，以换行符分隔
        cids = [int(cid.strip()) for cid in response.text.strip().split('\n') if cid.strip()]
        # 限制数量
        return cids[:limit] if len(cids) > limit else cids
    else:
        print(f"请求失败，状态码: {response.status_code}, 响应内容: {response.text}")
        return []


# 搜索含苯环的化合物
print("正在搜索含苯环的化合物...")
benzene_cids = search_benzene_compounds(20)
print(f"找到 {len(benzene_cids)} 个含苯环的化合物")
print(f"前10个CID: {benzene_cids[:10]}")


# 获取化合物详细信息
def get_compound_details(cids):
    """获取化合物的详细信息，包括名称和SMILES"""
    if not cids:
        return []
        
    # 将CID列表转换为字符串
    cid_str = ','.join(map(str, cids))
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_str}/property/Title,SMILES,InChI/JSON"
    
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        properties = data.get('PropertyTable', {}).get('Properties', [])
        return properties
    else:
        print(f"获取化合物详细信息失败，状态码: {response.status_code}")
        return []


# 获取化合物详细信息
print("正在获取含苯环化合物的详细信息...")
compound_details = get_compound_details(benzene_cids)
print(f"获取到 {len(compound_details)} 个化合物的详细信息")

# 创建DataFrame
df = pd.DataFrame(compound_details)
df['CID'] = benzene_cids
print(df.head())


# 初始化CC对象
print("正在初始化ComponentContribution对象...")
cc = ComponentContribution()
print("初始化完成")


# 计算pH=9下的生成自由能
print("开始计算pH=9下的生成自由能")
results = []

# 创建进度条
for i, compound in tqdm(enumerate(compound_details), total=len(compound_details), desc="计算进度"):
    try:
        # 使用InChI进行计算
        inchi = compound.get('InChI', '')
        if inchi:
            # 计算pH=9下的生成自由能
            dgf_prime, std = standard_dgf_prime_CC(inchi, cc, p_h=9.0, p_mg=3.0, I=0.25, T=298.15)
            results.append({
                'CID': benzene_cids[i],
                'Name': compound.get('Title', ''),
                'InChI': inchi,
                'SMILES': compound.get('SMILES', ''),
                'DGF_prime_pH9': dgf_prime.value.m,  # 获取数值部分
                'Uncertainty': std,
            })
        else:
            print(f"化合物 {compound.get('Title', '')} (CID: {benzene_cids[i]}) 没有InChI，跳过")
    except Exception as e:
        print(f"计算化合物 {compound.get('Title', '')} (CID: {benzene_cids[i]}) 时出错: {e}")
        # 添加空结果
        results.append({
            'CID': benzene_cids[i],
            'Name': compound.get('Title', ''),
            'InChI': '',
            'SMILES': compound.get('SMILES', ''),
            'DGF_prime_pH9': None,
            'Uncertainty': None,
        })

# 创建结果DataFrame
results_df = pd.DataFrame(results)
print(f"\n成功计算了 {len(results_df[results_df['DGF_prime_pH9'].notna()])} 个化合物的生成自由能")
print(results_df.head(10))

# 保存结果
output_path = 'benzene_ring_dgf_results_pH9.csv'
results_df.to_csv(output_path, index=False)
print(f"结果已保存到 {output_path}")

# 显示统计摘要
valid_results = results_df[results_df['DGF_prime_pH9'].notna()]
if not valid_results.empty:
    print(f"\n统计摘要:")
    print(f"成功计算的化合物数量: {len(valid_results)}")
    print(f"生成自由能范围: {valid_results['DGF_prime_pH9'].min():.2f} 到 {valid_results['DGF_prime_pH9'].max():.2f} kJ/mol")
    print(f"平均生成自由能: {valid_results['DGF_prime_pH9'].mean():.2f} kJ/mol")
    print(f"标准差: {valid_results['DGF_prime_pH9'].std():.2f} kJ/mol")