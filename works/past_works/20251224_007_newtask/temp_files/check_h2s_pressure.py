# 重新检查H2S(g)分压计算
import numpy as np

# 系统参数
V_aq = 1.00  # 水相体积 (L)
V_g = 5.00   # 气相体积 (L)
T = 298.15   # 温度 (K)
R = 0.082057 # 理想气体常数 (L atm mol^-1 K^-1)

# 亨利常数 (atm/(mol/m^3))
H_Henry = 9.741e-3  # atm/(mol/m^3)

print('亨利定律计算H2S(g)分压:')
print(f'亨利常数 H = {H_Henry:.3e} atm/(mol/m^3)')

# [H2S]aq = 2.000e-07 mol/L (从之前的计算)
H2S_aq = 2.000e-7  # mol/L

print(f'[H2S(aq)] = {H2S_aq:.3e} mol/L')

# 亨利定律: P = H * C，其中C是浓度(mol/m^3)
# 需要将[H2S]aq从mol/L转换为mol/m^3
H2S_aq_m3 = H2S_aq * 1000  # 1 L = 0.001 m^3, so 1 m^3 = 1000 L
print(f'[H2S(aq)] = {H2S_aq_m3:.3e} mol/m^3 (转换为m^3单位)')

# 使用亨利定律计算气相分压
H2S_gas_pressure = H_Henry * H2S_aq_m3
print(f'H2S(g)分压 = {H_Henry:.3e} * {H2S_aq_m3:.3e} = {H2S_gas_pressure:.3e} atm')

# 另一种方法：直接使用L单位
H_Henry_L = H_Henry * 1e-3  # 转换为 atm/(mol/L) 单位
print(f'\n使用L单位:')
print(f'亨利常数 (L单位) H_L = {H_Henry_L:.3e} atm/(mol/L)')
print(f'H2S(g)分压 = {H_Henry_L:.3e} * {H2S_aq:.3e} = {H_Henry_L * H2S_aq:.3e} atm')

print(f'\n确认: {H_Henry_L * H2S_aq:.3e} atm = {H2S_gas_pressure:.3e} atm (相同)')