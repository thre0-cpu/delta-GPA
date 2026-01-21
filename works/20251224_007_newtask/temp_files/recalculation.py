# 使用对数形式重新计算钼酸根水解平衡问题
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

print('平衡常数:')
print(f'K1 = {K1:.2e}, K2 = {K2:.2e}, K3 = {K3:.2e}, K4 = {K4:.2e}')
print(f'初始Mo总浓度 = {Mo_total:.2e} mol/L')

# 由于平衡常数很小，我们可以使用近似方法
# 假设水解程度很小，MoS4_2 ≈ Mo_total

# 设MoOS3_2 = x, 则根据平衡方程：
# K1 = [MoOS3_2][H2S]/[MoS4_2]
# K2 = [MoO2S2_2][H2S]/[MoOS3_2]
# K3 = [MoO3S_2][H2S]/[MoO2S2_2]
# K4 = [MoO4_2][H2S]/[MoO3S_2]

# 如果我们假设[H2S]aq ≈ constant，可以得到：
# [MoOS3_2] = K1 * [MoS4_2] / [H2S]
# [MoO2S2_2] = K2 * [MoOS3_2] / [H2S] = K1*K2 * [MoS4_2] / [H2S]^2
# [MoO3S_2] = K3 * [MoO2S2_2] / [H2S] = K1*K2*K3 * [MoS4_2] / [H2S]^3
# [MoO4_2] = K4 * [MoO3S_2] / [H2S] = K1*K2*K3*K4 * [MoS4_2] / [H2S]^4

# 同时，钼质量守恒：[MoS4_2] + [MoOS3_2] + [MoO2S2_2] + [MoO3S_2] + [MoO4_2] = Mo_total
# 气液平衡：[H2S]aq * V_aq + ([H2S]aq * H_Henry_L) * V_g = 总H2S摩尔数
# 其中总H2S摩尔数 = [MoOS3_2] + [MoO2S2_2] + [MoO3S_2] + [MoO4_2]

# 令 [H2S]aq = y，[MoS4_2] = z，则：
# [MoOS3_2] = K1 * z / y
# [MoO2S2_2] = K1*K2 * z / y^2
# [MoO3S_2] = K1*K2*K3 * z / y^3
# [MoO4_2] = K1*K2*K3*K4 * z / y^4

# 钼质量守恒: z + K1*z/y + K1*K2*z/y^2 + K1*K2*K3*z/y^3 + K1*K2*K3*K4*z/y^4 = Mo_total
# z * (1 + K1/y + K1*K2/y^2 + K1*K2*K3/y^3 + K1*K2*K3*K4/y^4) = Mo_total

# 气液平衡: y * V_aq + y * H_Henry_L * V_g = K1*z/y + K1*K2*z/y^2 + K1*K2*K3*z/y^3 + K1*K2*K3*K4*z/y^4
# y * (V_aq + H_Henry_L * V_g) = K1*z/y * (1 + K2/y + K2*K3/y^2 + K2*K3*K4/y^3)

# 这是一个复杂的非线性方程组，我们用迭代法求解

# 初始估计
y = 1e-10  # [H2S]aq的初始估计
z = Mo_total  # [MoS4_2]的初始估计

print(f'初始估计: [H2S]aq = {y:.2e}, [MoS4_2] = {z:.2e}')

for iteration in range(20):
    # 根据当前y和z计算其他浓度
    MoS4_2 = z
    MoOS3_2 = K1 * z / y
    MoO2S2_2 = K1*K2 * z / (y**2)
    MoO3S_2 = K1*K2*K3 * z / (y**3)
    MoO4_2 = K1*K2*K3*K4 * z / (y**4)
    
    # 检查钼质量守恒
    Mo_sum = MoS4_2 + MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
    
    # 根据气液平衡计算新的y值
    total_H2S_moles = MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
    y_new = total_H2S_moles / (V_aq + H_Henry_L * V_g)
    
    # 根据钼质量守恒计算新的z值
    z_new = Mo_total - (MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2)
    
    print(f'迭代 {iteration+1}: [H2S]aq = {y_new:.3e}, [MoS4_2] = {z_new:.3e}')
    
    # 检查收敛
    if abs(y_new - y) < 1e-20 and abs(z_new - z) < 1e-20:
        print('收敛完成')
        y = y_new
        z = z_new
        break
    
    y = y_new
    z = z_new

# 使用最终值计算所有浓度
MoS4_2 = z
MoOS3_2 = K1 * z / y
MoO2S2_2 = K1*K2 * z / (y**2)
MoO3S_2 = K1*K2*K3 * z / (y**3)
MoO4_2 = K1*K2*K3*K4 * z / (y**4)

# 计算H2S的其他存在形式
HS_neg = y * Ka1 / H_plus
S2_neg = HS_neg * Ka2 / H_plus

# 计算气相中H2S的分压
H2S_gas_conc = H_Henry_L * y
H2S_gas_pressure = H2S_gas_conc * R * T

print()
print('最终平衡浓度 (mol/L):')
print(f'[MoS4^2-] = {MoS4_2:.3e}')
print(f'[MoOS3^2-] = {MoOS3_2:.3e}')
print(f'[MoO2S2^2-] = {MoO2S2_2:.3e}')
print(f'[MoO3S^2-] = {MoO3S_2:.3e}')
print(f'[MoO4^2-] = {MoO4_2:.3e}')
print(f'[H2S(aq)] = {y:.3e}')
print(f'[HS^-] = {HS_neg:.3e}')
print(f'[S^2-] = {S2_neg:.3e}')
print(f'H2S(g)分压 = {H2S_gas_pressure:.3e} atm')

# 验证质量守恒
Mo_species_total = MoS4_2 + MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
print(f'\n钼质量守恒验证: {Mo_species_total:.3e} vs {Mo_total:.3e}')

total_H2S_moles = y * V_aq + H2S_gas_conc * V_g + HS_neg * V_aq + S2_neg * V_aq
expected_H2S_moles = MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
print(f'H2S质量守恒验证: {total_H2S_moles:.3e} vs {expected_H2S_moles:.3e}')

# 标准状态
c_std = 1.0  # mol/L
p_std = 1.0  # atm

# 含钼物种（浓度）
Mo_species = [MoS4_2, MoOS3_2, MoO2S2_2, MoO3S_2, MoO4_2]
Mo_species_names = ['MoS4^2-', 'MoOS3^2-', 'MoO2S2^2-', 'MoO3S^2-', 'MoO4^2-']

# 含硫不含钼物种（浓度或气压）
S_species = [y, HS_neg, S2_neg, H2S_gas_pressure]
S_species_names = ['H2S(aq)', 'HS^-(aq)', 'S^2-(aq)', 'H2S(g)']

# 计算归一化浓度
Mo_species_norm = [c/c_std for c in Mo_species]
S_species_norm = [c/c_std for c in S_species[:3]] + [H2S_gas_pressure/p_std]  # 最后一个是气压

print(f'\n归一化浓度:')
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