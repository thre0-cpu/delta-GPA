# 最终修正版 - 使用正确的平衡方法
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

# 由于平衡常数很小，水解程度应该很小
# 我们需要找到一个自洽的解，满足所有平衡条件和质量守恒

# 从平衡方程：
# [MoOS3_2] = K1 * [MoS4_2] / [H2S]
# [MoO2S2_2] = K2 * [MoOS3_2] / [H2S] = K1*K2 * [MoS4_2] / [H2S]^2
# [MoO3S_2] = K3 * [MoO2S2_2] / [H2S] = K1*K2*K3 * [MoS4_2] / [H2S]^3
# [MoO4_2] = K4 * [MoO3S_2] / [H2S] = K1*K2*K3*K4 * [MoS4_2] / [H2S]^4

# 钼质量守恒：
# [MoS4_2] + [MoOS3_2] + [MoO2S2_2] + [MoO3S_2] + [MoO4_2] = Mo_total

# 气液平衡：
# [H2S]aq * V_aq + ([H2S]aq * H_Henry_L) * V_g = [MoOS3_2] + [MoO2S2_2] + [MoO3S_2] + [MoO4_2]
# [H2S]aq * (V_aq + H_Henry_L * V_g) = [MoOS3_2] + [MoO2S2_2] + [MoO3S_2] + [MoO4_2]

# 令 [H2S]aq = x, [MoS4_2] = y，则：
# [MoOS3_2] = K1*y/x
# [MoO2S2_2] = K1*K2*y/x^2
# [MoO3S_2] = K1*K2*K3*y/x^3
# [MoO4_2] = K1*K2*K3*K4*y/x^4

# 钼质量守恒方程：
# y + K1*y/x + K1*K2*y/x^2 + K1*K2*K3*y/x^3 + K1*K2*K3*K4*y/x^4 = Mo_total
# y * (1 + K1/x + K1*K2/x^2 + K1*K2*K3/x^3 + K1*K2*K3*K4/x^4) = Mo_total
# y = Mo_total / (1 + K1/x + K1*K2/x^2 + K1*K2*K3/x^3 + K1*K2*K3*K4/x^4)

# 气液平衡方程：
# x * (V_aq + H_Henry_L * V_g) = K1*y/x + K1*K2*y/x^2 + K1*K2*K3*y/x^3 + K1*K2*K3*K4*y/x^4
# x * (V_aq + H_Henry_L * V_g) = y * (K1/x + K1*K2/x^2 + K1*K2*K3/x^3 + K1*K2*K3*K4/x^4)

# 将y的表达式代入：
# x * (V_aq + H_Henry_L * V_g) = [Mo_total / (1 + K1/x + K1*K2/x^2 + K1*K2*K3/x^3 + K1*K2*K3*K4/x^4)] * (K1/x + K1*K2/x^2 + K1*K2*K3/x^3 + K1*K2*K3*K4/x^4)

# 这是一个关于x的非线性方程，我们用数值方法求解
def equation(x):
    A = V_aq + H_Henry_L * V_g  # 1.000048705
    
    # 计算分母和分子
    denominator = 1 + K1/x + K1*K2/(x**2) + K1*K2*K3/(x**3) + K1*K2*K3*K4/(x**4)
    numerator = K1/x + K1*K2/(x**2) + K1*K2*K3/(x**3) + K1*K2*K3*K4/(x**4)
    
    left_side = x * A
    right_side = Mo_total * numerator / denominator
    
    return left_side - right_side

# 由于平衡常数很小，x应该也很小，我们尝试一些小值
import scipy.optimize as opt

try:
    # 使用brentq方法寻找根
    x_solution = opt.brentq(equation, 1e-15, 1e-5)
    print(f"找到 [H2S]aq = {x_solution:.3e} mol/L")
except:
    print("使用brentq未能找到解，尝试其他方法")
    # 尝试fsolve
    try:
        import scipy.optimize as opt
        x_solution = opt.fsolve(equation, 1e-10)[0]
        print(f"使用fsolve找到 [H2S]aq = {x_solution:.3e} mol/L")
    except:
        print("数值方法都失败，使用最简单近似")
        # 使用近似方法
        # 假设只有第一步反应显著发生
        x_solution = (K1 * Mo_total) ** 0.5
        print(f"近似 [H2S]aq = {x_solution:.3e} mol/L")

# 使用找到的[H2S]aq值计算其他浓度
x = x_solution

# 计算y = [MoS4_2]
denominator = 1 + K1/x + K1*K2/(x**2) + K1*K2*K3/(x**3) + K1*K2*K3*K4/(x**4)
y = Mo_total / denominator
MoS4_2 = y

# 计算其他钼物种
MoOS3_2 = K1 * y / x
MoO2S2_2 = K1*K2 * y / (x**2)
MoO3S_2 = K1*K2*K3 * y / (x**3)
MoO4_2 = K1*K2*K3*K4 * y / (x**4)

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