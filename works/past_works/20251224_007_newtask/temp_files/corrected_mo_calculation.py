# 重新计算钼酸根水解平衡问题 - 修正版
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

# 由于平衡常数很小，我们可以假设MoS4_2的浓度近似等于Mo_total
# 然后计算各步水解的产物浓度

# 设 [H2S]aq = x
# 从平衡方程可得：
# [MoOS3_2] = K1 * [MoS4_2] / [H2S] ≈ K1 * Mo_total / x
# [MoO2S2_2] = K2 * [MoOS3_2] / [H2S] ≈ K1*K2 * Mo_total / x^2
# [MoO3S_2] = K3 * [MoO2S2_2] / [H2S] ≈ K1*K2*K3 * Mo_total / x^3
# [MoO4_2] = K4 * [MoO3S_2] / [H2S] ≈ K1*K2*K3*K4 * Mo_total / x^4

# 由气液平衡：
# x * V_aq + (x * H_Henry_L) * V_g = [MoOS3_2] + [MoO2S2_2] + [MoO3S_2] + [MoO4_2]
# x * (V_aq + H_Henry_L * V_g) = K1*Mo_total/x + K1*K2*Mo_total/x^2 + K1*K2*K3*Mo_total/x^3 + K1*K2*K3*K4*Mo_total/x^4

# 令 A = V_aq + H_Henry_L * V_g
# A*x = K1*Mo_total/x + K1*K2*Mo_total/x^2 + K1*K2*K3*Mo_total/x^3 + K1*K2*K3*K4*Mo_total/x^4
# A*x^5 = K1*Mo_total*x^4 + K1*K2*Mo_total*x^3 + K1*K2*K3*Mo_total*x^2 + K1*K2*K3*K4*Mo_total*x
# A*x^4 = K1*Mo_total*x^3 + K1*K2*Mo_total*x^2 + K1*K2*K3*Mo_total*x + K1*K2*K3*K4*Mo_total

# 重新整理: A*x^4 - K1*Mo_total*x^3 - K1*K2*Mo_total*x^2 - K1*K2*K3*Mo_total*x - K1*K2*K3*K4*Mo_total = 0

# 由于K值很小，我们可以尝试用近似方法
# 将方程改写为: x^4 - (K1*Mo_total/A)*x^3 - (K1*K2*Mo_total/A)*x^2 - (K1*K2*K3*Mo_total/A)*x - (K1*K2*K3*K4*Mo_total/A) = 0

A = V_aq + H_Henry_L * V_g
print(f"A = V_aq + H_Henry_L * V_g = {A}")

# 计算多项式系数
coeffs = [1, -K1*Mo_total/A, -K1*K2*Mo_total/A, -K1*K2*K3*Mo_total/A, -K1*K2*K3*K4*Mo_total/A]
print(f"多项式系数: {coeffs}")

# 求解多项式方程
roots = np.roots(coeffs)
real_roots = [root.real for root in roots if np.isreal(root) and root.real > 0]

print(f"\n实数正根: {real_roots}")

if real_roots:
    x = min(real_roots)  # 选择最小的正实根
    print(f"取 [H2S]aq = {x:.3e} mol/L")
    
    # 现在我们使用正确的迭代方法：
    # 由于K值很小，我们可以用泰勒展开近似
    # 设MoS4_2 = Mo_total - a1, 其中a1是转化为MoOS3_2的量
    # 由于K很小，a1 << Mo_total，所以近似 [MoS4_2] ≈ Mo_total
    
    # 从平衡常数定义：
    # K1 = [MoOS3_2][H2S]/[MoS4_2] => [MoOS3_2] = K1*[MoS4_2]/[H2S]
    # K2 = [MoO2S2_2][H2S]/[MoOS3_2] => [MoO2S2_2] = K2*[MoOS3_2]/[H2S]
    # K3 = [MoO3S_2][H2S]/[MoO2S2_2] => [MoO3S_2] = K3*[MoO2S2_2]/[H2S]
    # K4 = [MoO4_2][H2S]/[MoO3S_2] => [MoO4_2] = K4*[MoO3S_2]/[H2S]
    
    # 使用近似 [MoS4_2] ≈ Mo_total
    MoS4_2 = Mo_total
    MoOS3_2 = K1 * MoS4_2 / x
    MoO2S2_2 = K2 * MoOS3_2 / x
    MoO3S_2 = K3 * MoO2S2_2 / x
    MoO4_2 = K4 * MoO3S_2 / x
    
    # 现在根据钼质量守恒修正MoS4_2
    total_other_Mo = MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
    if total_other_Mo < Mo_total:
        MoS4_2 = Mo_total - total_other_Mo
    else:
        # 如果其他钼物种总和超过了Mo_total，说明我们的[H2S]aq值过大
        # 我们需要重新求解
        print("需要重新求解，因为其他钼物种总和超过了Mo_total")
        # 使用迭代方法重新求解
        x_new = x
        for i in range(20):
            # 基于当前x_new计算钼物种
            MoOS3_2 = K1 * Mo_total / x_new  # 先用Mo_total近似
            MoO2S2_2 = K2 * MoOS3_2 / x_new
            MoO3S_2 = K3 * MoO2S2_2 / x_new
            MoO4_2 = K4 * MoO3S_2 / x_new
            
            # 根据钼质量守恒调整MoS4_2
            MoS4_2 = Mo_total - (MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2)
            if MoS4_2 < 0:
                MoS4_2 = 1e-30  # 防止负值
            
            # 根据气液平衡重新计算x_new
            total_H2S_moles = MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
            x_new_iter = total_H2S_moles / A
            
            print(f"迭代 {i+1}: x_new = {x_new_iter:.3e}, MoS4_2 = {MoS4_2:.3e}")
            
            # 检查收敛
            if abs(x_new_iter - x_new) < 1e-20:
                x_new = x_new_iter
                break
            x_new = x_new_iter
        
        # 使用最终值
        x = x_new
        MoOS3_2 = K1 * MoS4_2 / x
        MoO2S2_2 = K2 * MoOS3_2 / x
        MoO3S_2 = K3 * MoO2S2_2 / x
        MoO4_2 = K4 * MoO3S_2 / x
        
    print(f"\n最终钼物种浓度 (mol/L):")
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
        error = abs(Mo_species_total - Mo_total)/Mo_total
        print(f'钼平衡误差: {error:.2e}')
    else:
        print('钼平衡误差: 无法计算（Mo_total为0）')
    
    total_H2S_moles = x * V_aq + H2S_gas_conc * V_g + HS_neg * V_aq + S2_neg * V_aq
    expected_H2S_moles = MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
    print(f'H2S质量守恒验证: {total_H2S_moles:.3e} vs {expected_H2S_moles:.3e}')
    if expected_H2S_moles != 0:
        error = abs(total_H2S_moles - expected_H2S_moles)/expected_H2S_moles
        print(f'H2S平衡误差: {error:.2e}')
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
    if MoS4_2 != 0:
        calc_K1 = MoOS3_2 * x / MoS4_2
        if K1 != 0:
            error1 = abs(calc_K1-K1)/K1
            print(f'  实际K1 = {calc_K1:.3e}, 给定K1 = {K1:.3e}, 误差 = {error1:.2e}')
    else:
        print('  实际K1: 无法计算（MoS4_2为0）')
    
    if MoOS3_2 != 0:
        calc_K2 = MoO2S2_2 * x / MoOS3_2
        if K2 != 0:
            error2 = abs(calc_K2-K2)/K2
            print(f'  实际K2 = {calc_K2:.3e}, 给定K2 = {K2:.3e}, 误差 = {error2:.2e}')
    else:
        print('  实际K2: 无法计算（MoOS3_2为0）')
    
    if MoO2S2_2 != 0:
        calc_K3 = MoO3S_2 * x / MoO2S2_2
        if K3 != 0:
            error3 = abs(calc_K3-K3)/K3
            print(f'  实际K3 = {calc_K3:.3e}, 给定K3 = {K3:.3e}, 误差 = {error3:.2e}')
    else:
        print('  实际K3: 无法计算（MoO2S2_2为0）')
    
    if MoO3S_2 != 0:
        calc_K4 = MoO4_2 * x / MoO3S_2
        if K4 != 0:
            error4 = abs(calc_K4-K4)/K4
            print(f'  实际K4 = {calc_K4:.3e}, 给定K4 = {K4:.3e}, 误差 = {error4:.2e}')
    else:
        print('  实际K4: 无法计算（MoO3S_2为0）')
        
else:
    print("没有找到合适的正实根")