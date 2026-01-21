# 修正的测试equilibrator_api功能
from equilibrator_api import ComponentContribution, Q_, conditions
from equilibrator_cache import global_cache

print("开始测试equilibrator_api...")

# 创建ComponentContribution实例
print("正在创建ComponentContribution实例...")
cc = ComponentContribution()
print("ComponentContribution实例创建成功")

# 查看ComponentContribution对象的可用方法
print("ComponentContribution对象可用的方法和属性:")
methods = [method for method in dir(cc) if not method.startswith('_')]
print(methods)

# 重新尝试使用更合适的条件设置方法
try:
    # 使用条件对象设置
    my_ph = 9.0
    my_ionic_strength = Q_(0.25, 'M')
    my_temperature = Q_(298.15, 'K')
    
    # 创建条件
    c = conditions.PhysiologicalConditions(
        pH=my_ph,
        ionic_strength=my_ionic_strength,
        temperature=my_temperature
    )
    
    print(f"设置条件: pH={my_ph}, 离子强度={my_ionic_strength}, 温度={my_temperature}")
    
    # 定义反应
    reaction_formula = 'ATP + D-glucose = D-glucose 6-phosphate + ADP'
    print(f'正在解析反应: {reaction_formula}')

    # 解析反应
    parsed_rxn = cc.parse_reaction_formula(reaction_formula)
    print('反应解析完成')
    print('反应式:', str(parsed_rxn))
    
    # 使用指定条件计算标准反应自由能变化
    dg_prime = cc.standard_dg_prime(parsed_rxn, c)
    print(f'标准反应自由能变化: {dg_prime}')
    print(f'ΔG\'° 数值 (kJ/mol): {dg_prime.value.m:.2f}')
    
except Exception as e:
    print(f"使用PhysiologicalConditions仍然失败: {e}")
    
    # 尝试直接在计算时指定条件
    try:
        # 重新尝试解析反应
        reaction_formula = 'ATP + D-glucose = D-glucose 6-phosphate + ADP'
        print(f'正在解析反应: {reaction_formula}')
        parsed_rxn = cc.parse_reaction_formula(reaction_formula)
        print('反应解析完成')
        print('反应式:', str(parsed_rxn))
        
        # 使用cc的pH, ionic_strength, temperature参数
        from equilibrator_api import ComponentContribution
        from equilibrator_api.conditions import ConstantConcentrationCondition
        from equilibrator_api import Q_
        
        # 设置physiological条件
        phys_cond = cc.default_p_h, cc.default_ionic_strength, cc.default_temperature
        print(f"默认条件: pH={phys_cond[0]}, Ionic strength={phys_cond[1]}, Temperature={phys_cond[2]}")
        
        # 使用默认条件计算 (稍后更新条件)
        dg_prime = cc.standard_dg_prime(
            parsed_rxn,
            pH=9.0,
            ionic_strength=Q_(0.25, 'M'),
            temperature=Q_(298.15, 'K')
        )
        
        print(f'标准反应自由能变化: {dg_prime}')
        print(f'ΔG\'° 数值 (kJ/mol): {dg_prime.value.m:.2f}')
        
    except Exception as e2:
        print(f"第二种方法也失败: {e2}")
        print("尝试查看CC的README示例...")
        
        # 查看示例中的正确用法
        try:
            # 直接使用默认设置计算反应
            reaction_formula = 'ATP + D-glucose = D-glucose 6-phosphate + ADP'
            parsed_rxn = cc.parse_reaction_formula(reaction_formula)
            print('反应解析完成')
            print('反应式:', str(parsed_rxn))
            
            # 不指定特殊条件，使用默认条件计算
            dg_prime = cc.standard_dg_prime(
                parsed_rxn,
                pH=9.0,
                ionic_strength=Q_(0.25, 'M'),
                temperature=Q_(298.15, 'K')
            )
            
            print(f'标准反应自由能变化: {dg_prime}')
            print(f'ΔG\'° 数值 (kJ/mol): {dg_prime.value.m:.2f}')
            
        except Exception as e3:
            print(f"最终尝试也失败: {e3}")
            print("可能需要查阅CC_examples.ipynb中的具体示例")