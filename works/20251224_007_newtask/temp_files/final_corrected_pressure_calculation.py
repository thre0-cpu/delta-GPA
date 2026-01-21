# 最终修正版 - 修正H2S(g)分压计算（修复版）
import numpy as np

# 定义常数
K1 = 1.3e-5  # 第一个水解反应的平衡常数
K2 = 1.0e-5  # 第二个水解反应的平衡常数
K3 = 1.6e-5  # 第三个水解反应的平衡常数
K4 = 6.5e-6  # 第四个水解反应的平衡常数

# H2S的电离常数
pKa1 = 6.88
pKa2 = 12.90
Ka1 = 10**(-pKa1)  # = 1.318e-7
Ka2 = 10**(-pKa2)  # = 1.259e-13

# 系统参数
pH = 5.00
H_plus = 10**(-pH)  # [H+] = 1.0e-5 M
V_aq = 1.00  # 水相体积 (L)
V_g = 5.00   # 气相体积 (L)
T = 298.15   # 温度 (K)
R = 0.082057 # 理想气体常数 (L atm mol^-1 K^-1)

# 亨利常数 (atm/(mol/m^3))
H_Henry = 9.741e-3  # atm/(mol/m^3)

# 初始钼浓度
Mo_total = 2.00e-7  # mol/L

print('平衡常数和系统参数:')
print(f'K1 = {K1:.2e}, K2 = {K2:.2e}, K3 = {K3:.2e}, K4 = {K4:.2e}')
print(f'Ka1 = {Ka1:.2e}, Ka2 = {Ka2:.2e}')
print(f'pH = {pH}, [H+] = {H_plus:.2e} M')
print(f'V_aq = {V_aq} L, V_g = {V_g} L')
print(f'亨利常数 = {H_Henry:.2e} atm/(mol/m^3)')
print(f'初始Mo总浓度 = {Mo_total:.2e} mol/L')
print()

# 求解平衡
def equation(x):
    # A = V_aq + H_Henry_L * V_g，但需要正确处理单位
    # H_Henry是atm/(mol/m^3)，而浓度是mol/L
    # 需要将H_Henry转换为atm/(mol/L)：H_Henry_L = H_Henry / 1000
    H_Henry_L = H_Henry / 1000  # 因为1 m^3 = 1000 L
    
    A = V_aq + H_Henry_L * V_g
    
    # 计算分母和分子
    denominator = 1 + K1/x + K1*K2/(x**2) + K1*K2*K3/(x**3) + K1*K2*K3*K4/(x**4)
    numerator = K1/x + K1*K2/(x**2) + K1*K2*K3/(x**3) + K1*K2*K3*K4/(x**4)
    
    left_side = x * A
    right_side = Mo_total * numerator / denominator
    
    return left_side - right_side

# 使用数值方法求解
import scipy.optimize as opt

try:
    x_solution = opt.brentq(equation, 1e-15, 1e-3)
    print(f"找到 [H2S]aq = {x_solution:.3e} mol/L")
except:
    print("brentq方法失败，尝试fsolve")
    x_solution = opt.fsolve(equation, 1e-10)[0]
    print(f"使用fsolve找到 [H2S]aq = {x_solution:.3e} mol/L")

x = x_solution

# 计算钼物种浓度
H_Henry_L = H_Henry / 1000  # 从m^3转换为L单位
denominator = 1 + K1/x + K1*K2/(x**2) + K1*K2*K3/(x**3) + K1*K2*K3*K4/(x**4)
y = Mo_total / denominator
MoS4_2 = y

MoOS3_2 = K1 * y / x
MoO2S2_2 = K1*K2 * y / (x**2)
MoO3S_2 = K1*K2*K3 * y / (x**3)
MoO4_2 = K1*K2*K3*K4 * y / (x**4)

print(f"\n钼物种浓度 (mol/L):")
print(f'[MoS4^2-] = {MoS4_2:.3e}')
print(f'[MoOS3^2-] = {MoOS3_2:.3e}')
print(f'[MoO2S2^2-] = {MoO2S2_2:.3e}')
print(f'[MoO3S^2-] = {MoO3S_2:.3e}')
print(f'[MoO4^2-] = {MoO4_2:.3e}')

# 验证钼质量守恒
Mo_species_total = MoS4_2 + MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
print(f'\n钼质量守恒验证: {Mo_species_total:.3e} vs {Mo_total:.3e}')
if Mo_total != 0:
    error = abs(Mo_species_total - Mo_total)/Mo_total
    print(f'钼平衡误差: {error:.2e}')

# 计算含硫物种
HS_neg = Ka1 * x / H_plus
S2_neg = Ka2 * HS_neg / H_plus

print(f"\n含硫物种浓度 (mol/L):")
print(f'[H2S(aq)] = {x:.3e}')
print(f'[HS^-(aq)] = {HS_neg:.3e}')
print(f'[S^2-(aq)] = {S2_neg:.3e}')

# 修正H2S(g)分压计算
# 亨利定律: P = H * C，其中C是浓度(mol/m^3)
# [H2S]aq = x mol/L = x * 1000 mol/m^3
H2S_aq_m3 = x * 1000  # 转换为mol/m^3
H2S_gas_pressure = H_Henry * H2S_aq_m3
print(f'H2S(g)分压 (使用mol/m^3单位) = {H_Henry:.3e} * {H2S_aq_m3:.3e} = {H2S_gas_pressure:.3e} atm')

# 验证H2S质量守恒
# 水相中H2S摩尔数
H2S_aq_moles = x * V_aq
# 气相中H2S摩尔数 (通过理想气体定律 n = PV/RT)
H2S_gas_moles = H2S_gas_pressure * V_g / (R * T)
# 水相中HS-摩尔数
HS_moles = HS_neg * V_aq
# 水相中S2-摩尔数
S2_moles = S2_neg * V_aq

total_S_moles = H2S_aq_moles + H2S_gas_moles + HS_moles + S2_moles
expected_S_moles = MoOS3_2 + 2*MoO2S2_2 + 3*MoO3S_2 + 4*MoO4_2  # 每个钼物种含S数
expected_H2S_moles = MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2  # 仅考虑作为H2S产生的部分

print(f'\nH2S质量守恒验证:')
print(f'水相H2S: {H2S_aq_moles:.3e} mol')
print(f'气相H2S: {H2S_gas_moles:.3e} mol')
print(f'水相HS-: {HS_moles:.3e} mol')
print(f'水相S2-: {S2_moles:.3e} mol')
print(f'总S: {total_S_moles:.3e} mol')
print(f'由水解产生的H2S: {expected_H2S_moles:.3e} mol')

if expected_H2S_moles != 0:
    error = abs(H2S_aq_moles + H2S_gas_moles - expected_H2S_moles)/expected_H2S_moles
    print(f'H2S平衡误差: {error:.2e}')
else:
    print('H2S平衡误差: 无法计算（期望值为0）')

# 验证平衡常数
print(f'\n验证平衡常数:')
calc_K1 = MoOS3_2 * x / MoS4_2
calc_K2 = MoO2S2_2 * x / MoOS3_2
calc_K3 = MoO3S_2 * x / MoO2S2_2
calc_K4 = MoO4_2 * x / MoO3S_2

print(f'  实际K1 = {calc_K1:.3e}, 给定K1 = {K1:.3e}, 误差 = {abs(calc_K1-K1)/K1:.2e}')
print(f'  实际K2 = {calc_K2:.3e}, 给定K2 = {K2:.3e}, 误差 = {abs(calc_K2-K2)/K2:.2e}')
print(f'  实际K3 = {calc_K3:.3e}, 给定K3 = {K3:.3e}, 误差 = {abs(calc_K3-K3)/K3:.2e}')
print(f'  实际K4 = {calc_K4:.3e}, 给定K4 = {K4:.3e}, 误差 = {abs(calc_K4-K4)/K4:.2e}')

# 计算浓度/气压乘积比值
c_std = 1.0  # mol/L
p_std = 1.0  # atm

# 含钼物种（浓度）
Mo_species = [MoS4_2, MoOS3_2, MoO2S2_2, MoO3S_2, MoO4_2]
Mo_species_norm = [c/c_std for c in Mo_species]

# 含硫不含钼物种（浓度或气压）
S_species = [x, HS_neg, S2_neg, H2S_gas_pressure]
S_species_norm = [c/c_std for c in S_species[:3]] + [H2S_gas_pressure/p_std]

print(f'\n归一化浓度/气压:')
print(f'含钼物种: {[f"{val:.3e}" for val in Mo_species_norm]}')
print(f'含硫不含钼物种: {[f"{val:.3e}" for val in S_species_norm]}')

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