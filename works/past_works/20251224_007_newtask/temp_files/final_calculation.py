# 使用更精确的方法重新计算钼酸根水解平衡问题（修复版）
import numpy as np

# 定义常数
K1 = 1.3e-5  # 第一个水解反应的平衡常数
K2 = 1.0e-5  # 第二个水解反应的平衡常数
K3 = 1.6e-5  # 第三个水解反应的平衡常数
K4 = 6.5e-6  # 第四个水解反应的平衡常数

# H2S的电离常数
pKa1 = 6.88
pKa2 = 12.90
Ka1 = 10**(-pKa1)
Ka2 = 10**(-pKa2)

# 系统参数
pH = 5.00
H_plus = 10**(-pH)  # [H+] = 1.0e-5 M
V_aq = 1.00  # 水相体积 (L)
V_g = 5.00   # 气相体积 (L)
T = 298.15   # 温度 (K)
R = 0.082057 # 理想气体常数 (L atm mol^-1 K^-1)

# 亨利常数 (atm/(mol/m^3))
H_Henry = 9.741e-3  # atm/(mol/m^3)
# 转换为 atm/(mol/L) 单位
H_Henry_L = H_Henry * 1e-3  # 1 m^3 = 1000 L

# 初始钼浓度
Mo_total = 2.00e-7  # mol/L

print('平衡常数和系统参数:')
print(f'K1 = {K1:.2e}, K2 = {K2:.2e}, K3 = {K3:.2e}, K4 = {K4:.2e}')
print(f'初始Mo总浓度 = {Mo_total:.2e} mol/L')
print()

# 使用多项式方法求解
# A*x^4 - K1*Mo_total*x^3 - K1*K2*Mo_total*x^2 - K1*K2*K3*Mo_total*x - K1*K2*K3*K4*Mo_total = 0

# 计算系数
A = V_aq + H_Henry_L * V_g
coeffs = [A, -K1*Mo_total, -K1*K2*Mo_total, -K1*K2*K3*Mo_total, -K1*K2*K3*K4*Mo_total]

print(f"多项式系数: {coeffs}")

# 求解多项式方程
roots = np.roots(coeffs)
real_roots = [root.real for root in roots if np.isreal(root) and root.real > 0]

print(f"\n实数正根: {real_roots}")

if real_roots:
    x = min(real_roots)  # 选择最小的正实根
    print(f"取 [H2S]aq = {x:.3e} mol/L")
    
    # 计算各钼物种浓度
    MoS4_2 = Mo_total  # 近似认为未反应
    MoOS3_2 = K1 * MoS4_2 / x
    MoO2S2_2 = K2 * MoOS3_2 / x
    MoO3S_2 = K3 * MoO2S2_2 / x
    MoO4_2 = K4 * MoO3S_2 / x
    
    # 重新考虑钼质量守恒，调整MoS4_2浓度
    MoS4_2 = Mo_total - (MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2)
    if MoS4_2 < 0:
        print("\n警告：钼物种总浓度超过了初始浓度，需要重新调整计算方法")
        # 使用迭代方法
        MoS4_2 = Mo_total
        for i in range(10):
            MoOS3_2 = K1 * MoS4_2 / x
            MoO2S2_2 = K2 * MoOS3_2 / x
            MoO3S_2 = K3 * MoO2S2_2 / x
            MoO4_2 = K4 * MoO3S_2 / x
            MoS4_2_new = Mo_total - (MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2)
            if MoS4_2_new < 0:
                MoS4_2_new = 1e-30  # 防止负值
            if abs(MoS4_2_new - MoS4_2) < 1e-30:
                break
            MoS4_2 = MoS4_2_new
    
    print(f"\n调整后钼物种浓度 (mol/L):")
    print(f'[MoS4^2-] = {MoS4_2:.3e}')
    print(f'[MoOS3^2-] = {MoOS3_2:.3e}')
    print(f'[MoO2S2^2-] = {MoO2S2_2:.3e}')
    print(f'[MoO3S^2-] = {MoO3S_2:.3e}')
    print(f'[MoO4^2-] = {MoO4_2:.3e}')
    
    # 计算H2S的其他存在形式
    HS_neg = x * Ka1 / H_plus
    S2_neg = HS_neg * Ka2 / H_plus
    
    # 计算气相中H2S的分压
    H2S_gas_conc = H_Henry_L * x
    H2S_gas_pressure = H2S_gas_conc * R * T
    
    print(f"\n含硫物种浓度/气压 (mol/L 或 atm):")
    print(f'[H2S(aq)] = {x:.3e}')
    print(f'[HS^-] = {HS_neg:.3e}')
    print(f'[S^2-] = {S2_neg:.3e}')
    print(f'H2S(g)分压 = {H2S_gas_pressure:.3e}')
    
    # 验证质量守恒
    Mo_species_total = MoS4_2 + MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
    print(f'\n钼质量守恒验证: {Mo_species_total:.3e} vs {Mo_total:.3e}')
    if Mo_total != 0:
        print(f'钼平衡误差: {abs(Mo_species_total - Mo_total)/Mo_total:.2e}')
    else:
        print('钼平衡误差: 无法计算（Mo_total为0）')
    
    total_H2S_moles = x * V_aq + H2S_gas_conc * V_g + HS_neg * V_aq + S2_neg * V_aq
    expected_H2S_moles = MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
    print(f'H2S质量守恒验证: {total_H2S_moles:.3e} vs {expected_H2S_moles:.3e}')
    if expected_H2S_moles != 0:
        print(f'H2S平衡误差: {abs(total_H2S_moles - expected_H2S_moles)/expected_H2S_moles:.2e}')
    else:
        print('H2S平衡误差: 无法计算（期望值为0）')

    # 标准状态
    c_std = 1.0  # mol/L
    p_std = 1.0  # atm

    # 含钼物种（浓度）
    Mo_species = [MoS4_2, MoOS3_2, MoO2S2_2, MoO3S_2, MoO4_2]
    Mo_species_names = ['MoS4^2-', 'MoOS3^2-', 'MoO2S2^2-', 'MoO3S^2-', 'MoO4^2-']

    # 含硫不含钼物种（浓度或气压）
    S_species = [x, HS_neg, S2_neg, H2S_gas_pressure]
    S_species_names = ['H2S(aq)', 'HS^-(aq)', 'S^2-(aq)', 'H2S(g)']

    # 计算归一化浓度
    Mo_species_norm = [c/c_std for c in Mo_species]
    S_species_norm = [c/c_std for c in S_species[:3]] + [H2S_gas_pressure/p_std]  # 最后一个是气压

    print(f'\n归一化浓度/气压:')
    for name, val in zip(Mo_species_names, Mo_species_norm):
        print(f'{name}: {val:.3e}')
    for name, val in zip(S_species_names, S_species_norm):
        print(f'{name}: {val:.3e}')

    # 计算含钼物种浓度乘积
    Mo_product = 1.0
    for conc in Mo_species_norm:
        Mo_product *= conc

    # 计算含硫不含钼物种浓度/气压乘积
    S_product = 1.0
    for val in S_species_norm:
        S_product *= val

    # 计算比值
    ratio = Mo_product / S_product if S_product != 0 else float('inf')

    print(f'\n含钼物种浓度乘积: {Mo_product:.3e}')
    print(f'含硫不含钼物种浓度/气压乘积: {S_product:.3e}')
    print(f'浓度/气压乘积比值: {ratio:.3e}')
    
    # 验证平衡常数
    print(f'\n验证平衡常数:')
    calc_K1 = MoOS3_2 * x / MoS4_2 if MoS4_2 != 0 else float('inf')
    calc_K2 = MoO2S2_2 * x / MoOS3_2 if MoOS3_2 != 0 else float('inf')
    calc_K3 = MoO3S_2 * x / MoO2S2_2 if MoO2S2_2 != 0 else float('inf')
    calc_K4 = MoO4_2 * x / MoO3S_2 if MoO3S_2 != 0 else float('inf')

    print(f'  实际K1 = {calc_K1:.3e}, 给定K1 = {K1:.3e}, 误差 = {abs(calc_K1-K1)/K1:.2e if K1 != 0 else "N/A"}')
    print(f'  实际K2 = {calc_K2:.3e}, 给定K2 = {K2:.3e}, 误差 = {abs(calc_K2-K2)/K2:.2e if K2 != 0 else "N/A"}')
    print(f'  实际K3 = {calc_K3:.3e}, 给定K3 = {K3:.3e}, 误差 = {abs(calc_K3-K3)/K3:.2e if K3 != 0 else "N/A"}')
    print(f'  实际K4 = {calc_K4:.3e}, 给定K4 = {K4:.3e}, 误差 = {abs(calc_K4-K4)/K4:.2e if K4 != 0 else "N/A"}')
    
else:
    print("没有找到合适的正实根")