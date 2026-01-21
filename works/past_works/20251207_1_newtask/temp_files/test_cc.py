# 测试equilibrator_api功能
from equilibrator_api import ComponentContribution, Q_

print("开始测试equilibrator_api...")

# 创建ComponentContribution实例
print("正在创建ComponentContribution实例...")
cc = ComponentContribution()
print("ComponentContribution实例创建成功")

# 设置反应条件
temperature = Q_(298.15, 'K')
pH = 9.0
ionic_strength = Q_(0.25, 'M')

print(f"设置条件: pH={pH}, 离子强度={ionic_strength}, 温度={temperature}")

# 设置反应条件
cc.set_pH(pH)
cc.set_ionic_strength(ionic_strength)
cc.set_temperature(temperature)

print("条件设置完成")

# 定义反应
reaction_formula = 'ATP + 葡萄糖 = 葡萄糖-6-磷酸 + ADP'
print(f'正在解析反应: {reaction_formula}')

try:
    # 尝试使用英文名称，因为可能中文名无法识别
    reaction_formula_en = 'ATP + D-glucose = D-glucose 6-phosphate + ADP'
    parsed_rxn = cc.parse_reaction_formula(reaction_formula_en)
    print('反应解析完成')
    print('反应式:', str(parsed_rxn))
    
    # 计算标准反应自由能变化
    dg_prime = cc.standard_dg_prime(parsed_rxn)
    print(f'标准反应自由能变化: {dg_prime}')
    print(f'ΔG\'° 数值 (kJ/mol): {dg_prime.value.m:.2f}')
except Exception as e:
    print(f"解析反应时出错: {e}")
    print("尝试其他化合物名称格式...")
    
    # 尝试其他可能的名称格式
    try:
        reaction_formula_alt = 'aq:ATP + aq:D-glucose = aq:D-glucose-6-phosphate + aq:ADP'
        parsed_rxn = cc.parse_reaction_formula(reaction_formula_alt)
        print('使用替代格式解析成功')
        print('反应式:', str(parsed_rxn))
        
        # 计算标准反应自由能变化
        dg_prime = cc.standard_dg_prime(parsed_rxn)
        print(f'标准反应自由能变化: {dg_prime}')
        print(f'ΔG\'° 数值 (kJ/mol): {dg_prime.value.m:.2f}')
    except Exception as e2:
        print(f"使用替代格式仍然失败: {e2}")
        print("尝试KEGG ID格式...")
        
        # 尝试使用KEGG ID
        try:
            reaction_formula_kegg = 'C00002 + C00031 = C00092 + C00008'
            parsed_rxn = cc.parse_reaction_formula(reaction_formula_kegg)
            print('使用KEGG ID解析成功')
            print('反应式:', str(parsed_rxn))
            
            # 计算标准反应自由能变化
            dg_prime = cc.standard_dg_prime(parsed_rxn)
            print(f'标准反应自由能变化: {dg_prime}')
            print(f'ΔG\'° 数值 (kJ/mol): {dg_prime.value.m:.2f}')
        except Exception as e3:
            print(f"所有格式均失败: {e3}")
            print("请检查化合物名称是否正确，或者查看equilibrator_api文档以获取正确的命名方式")