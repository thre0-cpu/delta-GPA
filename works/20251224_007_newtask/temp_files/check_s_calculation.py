# 重新检查含硫物种计算
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
# 转换为 atm/(mol/L) 单位
H_Henry_L = H_Henry * 1e-3  # 1 m^3 = 1000 L

# 初始钼浓度
Mo_total = 2.00e-7  # mol/L

print('平衡常数和系统参数:')
print(f'K1 = {K1:.2e}, K2 = {K2:.2e}, K3 = {K3:.2e}, K4 = {K4:.2e}')
print(f'Ka1 = {Ka1:.2e}, Ka2 = {Ka2:.2e}')
print(f'pH = {pH}, [H+] = {H_plus:.2e} M')
print(f'初始Mo总浓度 = {Mo_total:.2e} mol/L')
print()

# 使用之前找到的解
x = 2.000e-07  # [H2S]aq

# 重新计算含硫物种
# H2S的电离平衡:
# H2S ⇌ H+ + HS-
# Ka1 = [H+][HS-]/[H2S]
# HS- ⇌ H+ + S2-
# Ka2 = [H+][S2-]/[HS-]

# [HS-] = Ka1 * [H2S] / [H+]
HS_neg = Ka1 * x / H_plus
print(f'使用公式 [HS-] = Ka1 * [H2S] / [H+]')
print(f'[HS-] = {Ka1:.2e} * {x:.2e} / {H_plus:.2e} = {HS_neg:.2e} mol/L')

# [S2-] = Ka2 * [HS-] / [H+]
S2_neg = Ka2 * HS_neg / H_plus
print(f'使用公式 [S2-] = Ka2 * [HS-] / [H+]')
print(f'[S2-] = {Ka2:.2e} * {HS_neg:.2e} / {H_plus:.2e} = {S2_neg:.2e} mol/L')

print(f'\nH2S电离分数:')
print(f'[H2S]: {x:.2e} mol/L')
print(f'[HS-]: {HS_neg:.2e} mol/L')
print(f'[S2-]: {S2_neg:.2e} mol/L')

# 检查电荷平衡近似
total_S_aq = x + HS_neg + S2_neg
print(f'\n水相中总硫浓度: {total_S_aq:.2e} mol/L')

# 计算电离度
alpha1 = HS_neg / (x + HS_neg + S2_neg)
alpha2 = S2_neg / (x + HS_neg + S2_neg)
print(f'H2S电离度 α1 = [HS-]/([H2S]+[HS-]+[S2-]) = {alpha1:.2e}')
print(f'H2S电离度 α2 = [S2-]/([H2S]+[HS-]+[S2-]) = {alpha2:.2e}')

# 计算气相中H2S的分压
H2S_gas_conc = H_Henry_L * x
H2S_gas_pressure = H2S_gas_conc * R * T

print(f'\n气相H2S分压 = {H2S_gas_pressure:.2e} atm')

# 计算含硫不含钼物种浓度/气压乘积
S_product = x * HS_neg * S2_neg * H2S_gas_pressure
print(f'\n含硫不含钼物种浓度/气压乘积 = [H2S(aq)] * [HS-(aq)] * [S2-(aq)] * P_H2S(g)')
print(f'= {x:.2e} * {HS_neg:.2e} * {S2_neg:.2e} * {H2S_gas_pressure:.2e}')
print(f'= {S_product:.2e}')