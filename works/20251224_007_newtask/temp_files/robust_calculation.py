# 使用更稳健的方法计算钼酸根水解平衡问题
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
print(f'Ka1 = {Ka1:.2e}, Ka2 = {Ka2:.2e}')
print(f'pH = {pH}, [H+] = {H_plus:.2e} M')
print(f'V_aq = {V_aq} L, V_g = {V_g} L')
print(f'亨利常数 = {H_Henry_L:.2e} atm/(mol/L)')
print(f'初始Mo总浓度 = {Mo_total:.2e} mol/L')

# 由于平衡常数很小，水解程度应该很小，我们使用近似方法
# 设x1, x2, x3, x4分别为各步反应的转化度
# [MoS4_2] = Mo_total - x1
# [MoOS3_2] = x1 - x2
# [MoO2S2_2] = x2 - x3
# [MoO3S_2] = x3 - x4
# [MoO4_2] = x4

# 平衡方程:
# K1 = [MoOS3_2][H2S]/[MoS4_2] = (x1-x2)*[H2S]/(Mo_total-x1)
# K2 = [MoO2S2_2][H2S]/[MoOS3_2] = (x2-x3)*[H2S]/(x1-x2)
# K3 = [MoO3S_2][H2S]/[MoO2S2_2] = (x3-x4)*[H2S]/(x2-x3)
# K4 = [MoO4_2][H2S]/[MoO3S_2] = x4*[H2S]/(x3-x4)

# 气液平衡: [H2S]aq * V_aq + ([H2S]aq * H_Henry_L) * V_g = x1 + x2 + x3 + x4
# [H2S]aq * (V_aq + H_Henry_L * V_g) = x1 + x2 + x3 + x4

# 由于x非常小，我们可以近似认为[MoS4_2] ≈ Mo_total
# 因此 K1 ≈ [MoOS3_2][H2S]/[MoS4_2] ≈ x1*[H2S]/Mo_total
# 所以 x1 ≈ K1*Mo_total/[H2S]

# 类似地，如果[H2S]已知，我们可以近似计算:
# x1 ≈ K1*Mo_total/[H2S]
# x2 ≈ K2*x1/[H2S] = K1*K2*Mo_total/[H2S]^2
# x3 ≈ K3*x2/[H2S] = K1*K2*K3*Mo_total/[H2S]^3
# x4 ≈ K4*x3/[H2S] = K1*K2*K3*K4*Mo_total/[H2S]^4

# 总H2S摩尔数 ≈ x1 + x2 + x3 + x4
# [H2S]aq * (V_aq + H_Henry_L * V_g) = x1 + x2 + x3 + x4

# 令 [H2S]aq = y，代入上面的公式：
# y * (V_aq + H_Henry_L * V_g) = K1*Mo_total/y + K1*K2*Mo_total/y^2 + K1*K2*K3*Mo_total/y^3 + K1*K2*K3*K4*Mo_total/y^4
# y * (V_aq + H_Henry_L * V_g) = K1*Mo_total * (1/y + K2/y^2 + K2*K3/y^3 + K2*K3*K4/y^4)

# 这是一个关于y的非线性方程，我们可以尝试求解
def equation(y):
    term1 = K1 * Mo_total / y
    term2 = K1 * K2 * Mo_total / (y**2)
    term3 = K1 * K2 * K3 * Mo_total / (y**3)
    term4 = K1 * K2 * K3 * K4 * Mo_total / (y**4)
    left_side = y * (V_aq + H_Henry_L * V_g)
    right_side = term1 + term2 + term3 + term4
    return left_side - right_side

# 由于平衡常数很小，[H2S]应该也很小，我们从一个小值开始尝试
y_values = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5]
results = [(y, abs(equation(y))) for y in y_values]
y_best = min(results, key=lambda x: x[1])[0]

print(f'\n寻找[H2S]aq的近似解:')
for y, err in results:
    print(f'  [H2S] = {y:.2e}, 误差 = {err:.2e}')
print(f'最佳近似解: [H2S]aq ≈ {y_best:.2e}')

# 使用找到的[H2S]aq值计算各物种浓度
y = y_best
x1 = K1 * Mo_total / y
x2 = K1 * K2 * Mo_total / (y**2)
x3 = K1 * K2 * K3 * Mo_total / (y**3)
x4 = K1 * K2 * K3 * K4 * Mo_total / (y**4)

# 计算各钼物种浓度
MoS4_2 = Mo_total - x1
MoOS3_2 = x1 - x2
MoO2S2_2 = x2 - x3
MoO3S_2 = x3 - x4
MoO4_2 = x4

print(f'\n各步转化度:')
print(f'x1 = {x1:.3e}')
print(f'x2 = {x2:.3e}')
print(f'x3 = {x3:.3e}')
print(f'x4 = {x4:.3e}')

print(f'\n钼物种浓度 (mol/L):')
print(f'[MoS4^2-] = {MoS4_2:.3e}')
print(f'[MoOS3^2-] = {MoOS3_2:.3e}')
print(f'[MoO2S2^2-] = {MoO2S2_2:.3e}')
print(f'[MoO3S^2-] = {MoO3S_2:.3e}')
print(f'[MoO4^2-] = {MoO4_2:.3e}')

# 计算H2S的其他存在形式
HS_neg = y * Ka1 / H_plus
S2_neg = HS_neg * Ka2 / H_plus

# 计算气相中H2S的分压
H2S_gas_conc = H_Henry_L * y
H2S_gas_pressure = H2S_gas_conc * R * T

print(f'\n含硫物种浓度/气压 (mol/L 或 atm):')
print(f'[H2S(aq)] = {y:.3e}')
print(f'[HS^-] = {HS_neg:.3e}')
print(f'[S^2-] = {S2_neg:.3e}')
print(f'H2S(g)分压 = {H2S_gas_pressure:.3e}')

# 验证质量守恒
Mo_species_total = MoS4_2 + MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
print(f'\n钼质量守恒验证: {Mo_species_total:.3e} vs {Mo_total:.3e}')
print(f'钼平衡误差: {abs(Mo_species_total - Mo_total)/Mo_total:.2e}')

total_H2S_moles = y * V_aq + H2S_gas_conc * V_g + HS_neg * V_aq + S2_neg * V_aq
expected_H2S_moles = x1 + x2 + x3 + x4
print(f'H2S质量守恒验证: {total_H2S_moles:.3e} vs {expected_H2S_moles:.3e}')
print(f'H2S平衡误差: {abs(total_H2S_moles - expected_H2S_moles)/expected_H2S_moles:.2e if expected_H2S_moles != 0 else 0}')

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
calc_K1 = MoOS3_2 * y / MoS4_2 if MoS4_2 != 0 else float('inf')
calc_K2 = MoO2S2_2 * y / MoOS3_2 if MoOS3_2 != 0 else float('inf')
calc_K3 = MoO3S_2 * y / MoO2S2_2 if MoO2S2_2 != 0 else float('inf')
calc_K4 = MoO4_2 * y / MoO3S_2 if MoO3S_2 != 0 else float('inf')

print(f'  实际K1 = {calc_K1:.3e}, 给定K1 = {K1:.3e}, 误差 = {abs(calc_K1-K1)/K1:.2e if K1 != 0 else 0}')
print(f'  实际K2 = {calc_K2:.3e}, 给定K2 = {K2:.3e}, 误差 = {abs(calc_K2-K2)/K2:.2e if K2 != 0 else 0}')
print(f'  实际K3 = {calc_K3:.3e}, 给定K3 = {K3:.3e}, 误差 = {abs(calc_K3-K3)/K3:.2e if K3 != 0 else 0}')
print(f'  实际K4 = {calc_K4:.3e}, 给定K4 = {K4:.3e}, 误差 = {abs(calc_K4-K4)/K4:.2e if K4 != 0 else 0}')