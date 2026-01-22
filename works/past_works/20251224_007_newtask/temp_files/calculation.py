# 让我们先简化问题，逐步计算
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

print('初始参数:')
print(f'Mo总浓度: {Mo_total:.2e} mol/L')
print(f'K1: {K1:.2e}, K2: {K2:.2e}, K3: {K3:.2e}, K4: {K4:.2e}')
print(f'pH: {pH}, [H+]: {H_plus:.2e} M')
print(f'V_aq: {V_aq} L, V_g: {V_g} L')
print(f'亨利常数: {H_Henry_L:.2e} atm/(mol/L)')

# 由于平衡常数较小，反应进行程度较低，我们从近似计算开始
# 假设大部分为MoS4_2，少量发生水解
MoS4_2 = Mo_total * 0.99  # 初始近似

# 根据平衡常数计算各步水解产物
H2S_aq = 1e-10  # 初始H2S浓度很小

# 迭代计算
for iteration in range(10):
    # 根据当前H2S浓度计算钼物种
    MoOS3_2 = K1 * MoS4_2 / (H2S_aq + 1e-15)  # 避免除零
    MoO2S2_2 = K2 * MoOS3_2 / (H2S_aq + 1e-15)
    MoO3S_2 = K3 * MoO2S2_2 / (H2S_aq + 1e-15)
    MoO4_2 = K4 * MoO3S_2 / (H2S_aq + 1e-15)
    
    # 检查钼质量守恒
    Mo_sum = MoS4_2 + MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
    
    # 计算H2S平衡浓度
    total_H2S_moles = MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2  # 由水解产生的H2S摩尔数
    # 气液平衡: H2S_aq * V_aq + (H2S_aq * H_Henry_L) * V_g = total_H2S_moles
    # H2S_aq * (V_aq + H_Henry_L * V_g) = total_H2S_moles
    H2S_aq_new = total_H2S_moles / (V_aq + H_Henry_L * V_g)
    
    print(f'迭代 {iteration+1}:')
    print(f'  MoS4_2: {MoS4_2:.2e}, MoOS3_2: {MoOS3_2:.2e}')
    print(f'  MoO2S2_2: {MoO2S2_2:.2e}, MoO3S_2: {MoO3S_2:.2e}, MoO4_2: {MoO4_2:.2e}')
    print(f'  Mo_sum: {Mo_sum:.2e}, Mo_total: {Mo_total:.2e}')
    print(f'  H2S_aq: {H2S_aq_new:.2e}')
    
    # 更新H2S浓度
    H2S_aq = H2S_aq_new
    
    # 根据钼质量守恒调整MoS4_2
    MoS4_2 = Mo_total - (MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2)
    if MoS4_2 < 0:
        MoS4_2 = 1e-20  # 确保非负
    
    # 检查收敛
    if abs(H2S_aq_new - H2S_aq) < 1e-15 and abs(Mo_total - Mo_sum) < 1e-15:
        print('收敛完成')
        break

# 计算H2S的其他存在形式
HS_neg = H2S_aq * Ka1 / H_plus
S2_neg = HS_neg * Ka2 / H_plus

# 计算气相中H2S的分压
H2S_gas_conc = H_Henry_L * H2S_aq
H2S_gas_pressure = H2S_gas_conc * R * T

print('\n最终平衡浓度 (mol/L):')
print(f'[MoS4^2-] = {MoS4_2:.2e}')
print(f'[MoOS3^2-] = {MoOS3_2:.2e}')
print(f'[MoO2S2^2-] = {MoO2S2_2:.2e}')
print(f'[MoO3S^2-] = {MoO3S_2:.2e}')
print(f'[MoO4^2-] = {MoO4_2:.2e}')
print(f'[H2S(aq)] = {H2S_aq:.2e}')
print(f'[HS^-] = {HS_neg:.2e}')
print(f'[S^2-] = {S2_neg:.2e}')
print(f'H2S(g)分压 = {H2S_gas_pressure:.2e} atm')

# 验证质量守恒
Mo_species_total = MoS4_2 + MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
print(f'\n钼质量守恒验证: {Mo_species_total:.2e} vs {Mo_total:.2e}')

total_H2S_moles = H2S_aq * V_aq + H2S_gas_conc * V_g + HS_neg * V_aq + S2_neg * V_aq
expected_H2S_moles = MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
print(f'H2S质量守恒验证: {total_H2S_moles:.2e} vs {expected_H2S_moles:.2e}')

# 标准状态
c_std = 1.0  # mol/L
p_std = 1.0  # atm

# 含钼物种（浓度）
Mo_species = [MoS4_2, MoOS3_2, MoO2S2_2, MoO3S_2, MoO4_2]
Mo_species_names = ['MoS4^2-', 'MoOS3^2-', 'MoO2S2^2-', 'MoO3S^2-', 'MoO4^2-']

# 含硫不含钼物种（浓度或气压）
S_species = [H2S_aq, HS_neg, S2_neg, H2S_gas_pressure]
S_species_names = ['H2S(aq)', 'HS^-(aq)', 'S^2-(aq)', 'H2S(g)']

# 计算归一化浓度
Mo_species_norm = [c/c_std for c in Mo_species]
S_species_norm = [c/c_std for c in S_species[:3]] + [H2S_gas_pressure/p_std]  # 最后一个是气压

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

print('\n归一化浓度/气压:')
for name, val in zip(Mo_species_names, Mo_species_norm):
    print(f'{name}: {val:.2e}')
for name, val in zip(S_species_names, S_species_norm):
    print(f'{name}: {val:.2e}')

print(f'\n含钼物种浓度乘积: {Mo_product:.2e}')
print(f'含硫不含钼物种浓度/气压乘积: {S_product:.2e}')
print(f'浓度/气压乘积比值: {ratio:.2e}')