# 导入必要的库
import numpy as np
import pandas as pd
from equilibrator_api import ComponentContribution, Q_
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# 抑制RDKit警告
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

print("正在初始化ComponentContribution...")
cc = ComponentContribution()
print("ComponentContribution初始化完成")

# 尝试查找UDP化合物
compound_name = "UDP"

try:
    # 尝试通过不同的方式查找UDP
    udp_compound = cc.get_compound(compound_name)
    if udp_compound is None:
        udp_compound = cc.get_compound(f"kegg:{compound_name}")
    if udp_compound is None:
        udp_compound = cc.get_compound(f"chebi:{compound_name}")
    if udp_compound is None:
        udp_compound = cc.get_compound_by_inchi(compound_name)
    
    if udp_compound is not None:
        print(f"成功找到UDP化合物: {udp_compound}")
        print(f"ID: {udp_compound.id}")
        print(f"InChI: {udp_compound.inchi}")
        print(f"InChI Key: {udp_compound.inchi_key}")
        print(f"名称: {udp_compound.name}")
    else:
        print("未能找到UDP化合物，请尝试其他名称或格式")
        # 尝试其他可能的UDP别名
        possible_names = ["uridine_diphosphate", "uridine_5'-diphosphate", "UMP", "UTP", "ADP", "GDP"]
        for name in possible_names:
            compound = cc.get_compound(name)
            if compound is not None:
                print(f"找到类似化合物 {name}: {compound}")
except Exception as e:
    print(f"查找UDP化合物时出错: {e}")

# 如果上面没找到UDP，尝试更具体的搜索
if udp_compound is None:
    # 尝试使用KEGG ID
    try:
        udp_compound = cc.get_compound("kegg:C00015")  # UDP的KEGG ID
        if udp_compound is not None:
            print(f"通过KEGG ID找到UDP: {udp_compound.name}")
    except:
        pass

    # 如果还是没找到，尝试其他方式
    if udp_compound is None:
        try:
            udp_compound = cc.get_compound("uridine diphosphate")
            if udp_compound is not None:
                print(f"通过全名找到UDP: {udp_compound.name}")
        except:
            pass

if udp_compound is None:
    print("仍然无法找到UDP化合物，尝试使用更广泛的搜索")
else:
    print(f"成功获取UDP化合物信息: {udp_compound}")

# 计算UDP在指定条件下的生成自由能
if udp_compound is not None:
    # 设定条件
    temp = 298  # K
    ph = 10
    pmg = 6
    ionic_strength = 0.3  # M
    
    print(f"\n开始计算UDP在以下条件下的生成自由能:")
    print(f"温度: {temp} K")
    print(f"pH: {ph}")
    print(f"pMg: {pmg}")
    print(f"离子强度: {ionic_strength} M")
    
    try:
        # 使用transform方法将标准条件转换到指定条件
        transformed_compound = udp_compound.transform(p_h=Q_(ph), p_mg=Q_(pmg), ionic_strength=Q_(f'{ionic_strength}M'), temperature=Q_(f'{temp}K'))
        
        # 获取标准生成自由能
        standard_dgf = udp_compound.standard_dgf
        
        # 计算在指定条件下的生成自由能
        adjusted_dgf = standard_dgf + transformed_compound
        
        print(f"\nUDP在指定条件下的生成自由能 (ΔfG'°): {adjusted_dgf.value} kJ/mol")
        print(f"标准生成自由能 (ΔfG°): {standard_dgf.value} kJ/mol")
        print(f"校正项: {transformed_compound.value} kJ/mol")
        
        # 获取不确定性
        try:
            uncertainty = udp_compound.uncertainty
            print(f"不确定性: ±{uncertainty} kJ/mol")
        except:
            print("未找到不确定性数据")
            
    except Exception as e:
        print(f"计算过程中出现错误: {e}")
        
        # 尝试另一种方法
        try:
            print("\n尝试使用另一种方法计算...")
            # 使用compound的transformed formation energy方法
            formation_energy = udp_compound.transformed_formation_energy(p_h=Q_(ph), p_mg=Q_(pmg), ionic_strength=Q_(f'{ionic_strength}M'), temperature=Q_(f'{temp}K'))
            print(f"UDP在指定条件下的生成自由能 (ΔfG'°): {formation_energy.value} kJ/mol")
        except Exception as e2:
            print(f"第二种方法也失败了: {e2}")
else:
    print("由于找不到UDP化合物，无法进行计算")