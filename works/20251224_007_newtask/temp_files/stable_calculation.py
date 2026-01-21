# 使用最简单和最稳定的方法重新计算钼酸根水解平衡问题
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

# 由于平衡常数很小，我们可以假设水解程度非常小
# 因此 [MoS4^2-] ≈ Mo_total
# 我们可以将问题简化为一个简单的近似计算

# 假设 [H2S]aq = x，这是一个待求解的值
# 从平衡关系可得（近似[MoS4_2] ≈ Mo_total）：
# [MoOS3_2] ≈ K1 * Mo_total / x
# [MoO2S2_2] ≈ K1*K2 * Mo_total / x^2
# [MoO3S_2] ≈ K1*K2*K3 * Mo_total / x^3
# [MoO4_2] ≈ K1*K2*K3*K4 * Mo_total / x^4

# 由气液平衡：
# x * V_aq + (x * H_Henry_L) * V_g = [MoOS3_2] + [MoO2S2_2] + [MoO3S_2] + [MoO4_2]
# x * (V_aq + H_Henry_L * V_g) = K1*Mo_total/x + K1*K2*Mo_total/x^2 + K1*K2*K3*Mo_total/x^3 + K1*K2*K3*K4*Mo_total/x^4

# 令 A = V_aq + H_Henry_L * V_g
# A*x = K1*Mo_total/x + K1*K2*Mo_total/x^2 + K1*K2*K3*Mo_total/x^3 + K1*K2*K3*K4*Mo_total/x^4
# A*x^2 = K1*Mo_total + K1*K2*Mo_total/x + K1*K2*K3*Mo_total/x^2 + K1*K2*K3*K4*Mo_total/x^3

# 由于平衡常数非常小，[H2S]aq 也会很小，但不会为0
# 我们可以使用另一种方法：假设水解程度很小，可以近似计算

# 首先，由于K值很小，我们预期x（[H2S]aq）也会很小
# 从多项式求解得到的x值是 2.301e-07 mol/L

x = 2.301e-07  # [H2S]aq的近似值

print(f'使用近似 [H2S]aq = {x:.3e} mol/L')

# 计算各钼物种浓度，考虑钼质量守恒
# 使用迭代方法，但更稳定
MoS4_2 = Mo_total  # 初始近似

for i in range(5):  # 进行几次迭代
    MoOS3_2 = K1 * MoS4_2 / x
    MoO2S2_2 = K2 * MoOS3_2 / x
    MoO3S_2 = K3 * MoO2S2_2 / x
    MoO4_2 = K4 * MoO3S_2 / x
    
    total_other_Mo = MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
    
    if total_other_Mo < Mo_total:
        MoS4_2 = Mo_total - total_other_Mo
    else:
        # 如果其他钼物种总和超过了Mo_total，需要调整x
        # 这意味着我们对[H2S]aq的初始估计过大
        # 重新计算x值，基于Mo_total
        total_H2S_moles = total_other_Mo  # 由水解产生的H2S摩尔数
        x = total_H2S_moles / (V_aq + H_Henry_L * V_g)
        MoS4_2 = Mo_total  # 重新开始
        continue  # 重新迭代
    
    print(f'迭代 {i+1}: [MoS4_2] = {MoS4_2:.3e}, [H2S]aq = {x:.3e}')

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