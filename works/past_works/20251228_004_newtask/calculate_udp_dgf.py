from equilibrator_api import ComponentContribution, Q_, Reaction
import numpy as np

print("开始初始化ComponentContribution对象...")
cc = ComponentContribution()
print("ComponentContribution对象初始化完成")

# 设置条件
p_h = 10.0
p_mg = 6.0
T = 298.15
I = 0.25

# 设置CC对象的参数
cc.p_h = Q_(p_h)
cc.p_mg = Q_(p_mg)
cc.ionic_strength = Q_(f'{I}M')
cc.temperature = Q_(f'{T}K')

print(f"已设置条件: pH={p_h}, pMg={p_mg}, T={T}K, I={I}M")

# 尝试获取UDP化合物
udp_identifiers = ['UDP', 'kegg:C00015', 'uridine diphosphate', 'Uridine 5\'-diphosphate']

for identifier in udp_identifiers:
    print(f'尝试使用标识符: {identifier}')
    try:
        # 尝试获取化合物
        compound = None
        
        # 尝试不同的获取化合物方法
        methods = [
            lambda x: cc.get_compound(x),
            lambda x: cc.search_compound(x),
            lambda x: cc.get_compound(f'kegg:{x}'),
            lambda x: cc.get_compound_by_inchi(x) if x.startswith('InChI=') else None
        ]
        
        for method in methods:
            try:
                result = method(identifier)
                if result is not None:
                    compound = result
                    break
            except:
                continue
        
        if compound is not None:
            print(f'成功找到化合物: {compound}')
            
            # 创建虚拟反应: 0 -> 1 化合物
            rxn_c = Reaction({compound: 1})
            
            # 计算能量
            print("正在计算生成自由能...")
            dg_prime_measurement = cc.standard_dg_prime(rxn_c)
            
            # 提取数值
            val = dg_prime_measurement.value.m_as('kJ/mol')
            err = dg_prime_measurement.error.m_as('kJ/mol')
            
            print(f'UDP Δ_fG\'° 在 {T}K, pH={p_h}, pMg={p_mg} 条件下为: {val:.2f} ± {err:.2f} kJ/mol')
            break
        else:
            print('  无法找到该化合物')
    except Exception as e:
        print(f'  计算出错: {str(e)}')
else:
    print('无法找到UDP或计算其生成自由能')