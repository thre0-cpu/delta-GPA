# 使用正确方法重新计算钼酸根与硫代钼酸根水解平衡问题
import numpy as np
import scipy.optimize as opt

# 定义常数
K1 = 1.3e-5  # 第一个水解反应的平衡常数
K2 = 1.0e-5  # 第二个水解反应的平衡常数
K3 = 1.6e-5  # 第三个水解反应的平衡常数
K4 = 6.5e-6  # 第四个水解反应的平衡常数

# H2S的电离常数
pKa1 = 6.88
pKa2 = 12.90

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
H_Henry_L = H_Henry / 1000  # 因为1 m^3 = 1000 L

# 初始钼浓度
C0 = 2.00e-7  # mol/L

print('平衡常数和系统参数:')
print(f'K1 = {K1:.2e}, K2 = {K2:.2e}, K3 = {K3:.2e}, K4 = {K4:.2e}')
print(f'pH = {pH}, [H+] = {H_plus:.2e} M')
print(f'V_aq = {V_aq} L, V_g = {V_g} L')
print(f'亨利常数 = {H_Henry_L:.2e} atm/(mol/L)')
print(f'初始Mo总浓度 = {C0:.2e} mol/L')
print()

# 正确的计算方法
# 定义 α1 和 α2
alpha1 = 10**(pH - pKa1)  # = 10^(5.00-6.88) = 10^(-1.88) = 0.01318
alpha2 = 10**(2*pH - pKa1 - pKa2)  # = 10^(10.00-6.88-12.90) = 10^(-9.78) = 1.66e-10

print(f'α1 = {alpha1:.3e}, α2 = {alpha2:.3e}')

# 计算 F
F = (1 + alpha1 + alpha2) + H_Henry_L * V_g / (R * T)
print(f'(1 + α1 + α2) = {1 + alpha1 + alpha2:.3f}')
print(f'H_L * V_g / (R * T) = {H_Henry_L * V_g / (R * T):.3f}')
print(f'F = {F:.3f}')

# 定义 D(x) 和 S(x)
def D(x):
    return 1 + K1/x + K1*K2/(x**2) + K1*K2*K3/(x**3) + K1*K2*K3*K4/(x**4)

def S(x):
    return K1/x + 2*K1*K2/(x**2) + 3*K1*K2*K3/(x**3) + 4*K1*K2*K3*K4/(x**4)

# 求解方程 x = C0 * S(x) / (D(x) * F)
def equation(x):
    return x - C0 * S(x) / (D(x) * F)

# 使用数值方法求解
try:
    x_solution = opt.brentq(equation, 1e-10, 1e-6)
    print(f"找到 [H2S]aq = {x_solution:.4e} mol/L")
except:
    print("brentq方法失败，尝试fsolve")
    x_solution = opt.fsolve(equation, 1e-7)[0]
    print(f"使用fsolve找到 [H2S]aq = {x_solution:.4e} mol/L")

x = x_solution

# 计算钼物种浓度
A0 = C0 / D(x)  # [MoS4^2-]
A1 = K1 * A0 / x  # [MoOS3^2-]
A2 = K1*K2 * A0 / (x**2)  # [MoO2S2^2-]
A3 = K1*K2*K3 * A0 / (x**3)  # [MoO3S^2-]
A4 = K1*K2*K3*K4 * A0 / (x**4)  # [MoO4^2-]

print(f"\n钼物种浓度 (mol/L):")
print(f'[MoS4^2-] = {A0:.4e}')
print(f'[MoOS3^2-] = {A1:.4e}')
print(f'[MoO2S2^2-] = {A2:.4e}')
print(f'[MoO3S^2-] = {A3:.4e}')
print(f'[MoO4^2-] = {A4:.4e}')

# 验证钼质量守恒
total_Mo = A0 + A1 + A2 + A3 + A4
print(f'\n钼质量守恒验证: {total_Mo:.4e} vs {C0:.4e}')
if C0 != 0:
    error = abs(total_Mo - C0)/C0
    print(f'钼平衡误差: {error:.2e}')

# 计算含硫物种
HS_neg = alpha1 * x
S2_neg = alpha2 * x

H2S_gas_pressure = H_Henry * x * 1000  # 从mol/m^3转换
print(f"\n含硫物种浓度 (mol/L):")
print(f'[H2S(aq)] = {x:.4e}')
print(f'[HS^-(aq)] = {HS_neg:.4e}')
print(f'[S^2-(aq)] = {S2_neg:.4e}')
print(f'H2S(g)分压 = {H2S_gas_pressure:.4e} atm')

# 验证平衡常数
print(f'\n验证平衡常数:')
calc_K1 = A1 * x / A0
calc_K2 = A2 * x / A1
calc_K3 = A3 * x / A2
calc_K4 = A4 * x / A3

print(f'  实际K1 = {calc_K1:.3e}, 给定K1 = {K1:.3e}, 误差 = {abs(calc_K1-K1)/K1:.2e}')
print(f'  实际K2 = {calc_K2:.3e}, 给定K2 = {K2:.3e}, 误差 = {abs(calc_K2-K2)/K2:.2e}')
print(f'  实际K3 = {calc_K3:.3e}, 给定K3 = {K3:.3e}, 误差 = {abs(calc_K3-K3)/K3:.2e}')
print(f'  实际K4 = {calc_K4:.3e}, 给定K4 = {K4:.3e}, 误差 = {abs(calc_K4-K4)/K4:.2e}')

# 计算浓度/气压乘积比值
c_std = 1.0  # mol/L
p_std = 1.0  # atm

# 含钼物种（浓度）
Mo_species = [A0, A1, A2, A3, A4]
Mo_species_norm = [c/c_std for c in Mo_species]

# 含硫不含钼物种（浓度或气压）
S_species = [x, HS_neg, S2_neg, H2S_gas_pressure]
S_species_norm = [c/c_std for c in S_species[:3]] + [H2S_gas_pressure/p_std]

print(f'\n归一化浓度/气压:')
Mo_names = ['MoS4^2-', 'MoOS3^2-', 'MoO2S2^2-', 'MoO3S^2-', 'MoO4^2-']
S_names = ['H2S(aq)', 'HS^-(aq)', 'S^2-(aq)', 'H2S(g)']
for name, val in zip(Mo_names, Mo_species_norm):
    print(f'{name}: {val:.3e}')
for name, val in zip(S_names, S_species_norm):
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