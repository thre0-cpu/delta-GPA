from equilibrator_api import ComponentContribution, Q_
import numpy as np

print('开始初始化ComponentContribution...')
cc = ComponentContribution()
print('初始化完成')

# 设置条件
p_h = 10
p_mg = 6
T = 298  # K
I = 0.25  # M

cc.p_h = Q_(p_h)
cc.p_mg = Q_(p_mg)
cc.ionic_strength = Q_(f'{I}M')
cc.temperature = Q_(f'{T}K')

print(f'已设置条件: pH={p_h}, pMg={p_mg}, T={T}K, I={I}M')

# 尝试直接搜索UDP
print('正在搜索UDP化合物...')
udp = cc.search_compound('UDP')
if udp:
    print(f'找到UDP: {udp.id} (InChI: {udp.inchi})')
    print('正在计算生成自由能...')
    try:
        dgf, sigma_fin, sigma_inf = cc.standard_dg_formation(udp)
        std = np.linalg.norm(sigma_fin) if sigma_fin is not None else 0
        print(f'UDP在指定条件下的生成自由能为: {dgf:.2f} kJ/mol ± {std:.2f} kJ/mol')
    except Exception as e:
        print(f'计算生成自由能时出现错误: {e}')
        import traceback
        traceback.print_exc()
else:
    print('未能找到UDP化合物')
    # 尝试其他可能的标识符
    print('尝试使用KEGG ID C00015...')
    try:
        udp = cc.get_compound('kegg:C00015')
        if udp:
            print(f'找到UDP: {udp.id} (InChI: {udp.inchi})')
            print('正在计算生成自由能...')
            dgf, sigma_fin, sigma_inf = cc.standard_dg_formation(udp)
            std = np.linalg.norm(sigma_fin) if sigma_fin is not None else 0
            print(f'UDP在指定条件下的生成自由能为: {dgf:.2f} kJ/mol ± {std:.2f} kJ/mol')
        else:
            print('使用KEGG ID也未能找到UDP')
    except Exception as e:
        print(f'尝试KEGG ID时出现错误: {e}')