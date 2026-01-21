# ATP + 葡萄糖反应热力学平衡计算脚本
import subprocess
import json
import sys
import os
from typing import Optional, Tuple, Union

import numpy as np
import numpy.typing as npt
from rdkit import Chem

from equilibrator_api import ComponentContribution, Q_
from scipy.optimize import fsolve
from scipy.constants import R
import math

print("正在初始化ComponentContribution实例...")

# 初始化CC类
cc = ComponentContribution()
print("ComponentContribution实例已创建")

# 设置反应条件
p_h = 9.0
p_mg = 3.0  # 默认值
I = 0.25  # 离子强度 (M)
T = 298.15  # 温度 (K)

cc.p_h = Q_(p_h)
cc.p_mg = Q_(p_mg)
cc.ionic_strength = Q_(f'{I}M')
cc.temperature = Q_(f'{T}K')

print(f'设置条件: pH={p_h}, pMg={p_mg}, 离子强度={I}M, 温度={T}K')

# 定义反应，使用KEGG ID
# ATP: C00002, D-glucose: C00031, D-glucose 6-phosphate: C00092, ADP: C00008
reaction_formula = 'C00002 + C00031 = C00092 + C00008'
print(f'正在解析反应: {reaction_formula}')

# 解析反应
parsed_rxn = cc.parse_reaction_formula(reaction_formula)
print('反应解析完成')
print('反应式:', str(parsed_rxn))

# 计算标准反应自由能变化
dg_prime = cc.standard_dg_prime(parsed_rxn)
print(f'标准反应自由能变化: {dg_prime}')
print(f'ΔG\'° 数值 (kJ/mol): {dg_prime.value.m:.2f}')

# 从反应自由能计算平衡常数
T = 298.15  # K
delta_g_prime_kj_per_mol = dg_prime.value.m  # kJ/mol
delta_g_prime_j_per_mol = delta_g_prime_kj_per_mol * 1000  # J/mol

# 计算平衡常数 K_eq
RT = R * T
ln_keq = -delta_g_prime_j_per_mol / RT
K_eq = math.exp(ln_keq)

print(f"平衡常数 K_eq = {K_eq:.2e}")

# 计算平衡时各组分浓度
# 对于反应 ATP + 葡萄糖 = 葡萄糖-6-磷酸 + ADP
# 设反应进行的程度为 x
# [ATP] = 0.001 - x
# [葡萄糖] = 0.001 - x
# [葡萄糖-6-磷酸] = x
# [ADP] = x

# 平衡方程: K_eq = ([G6P][ADP])/([ATP][glucose]) = x^2 / ((0.001-x)^2)

# 定义方程求解x
def equilibrium_equation(x):
    numerator = x**2
    denominator = (0.001 - x)**2
    
    # 确保不除以0
    if denominator == 0:
        return float('inf') if numerator > 0 else 0
    
    return numerator / denominator - K_eq

# 初始猜测值
initial_guess = 0.0001  # 小于初始浓度0.001

# 求解反应进行程度x
x_solution = fsolve(equilibrium_equation, initial_guess)[0]

# 确保x在合理范围内
if x_solution < 0 or x_solution > 0.001:
    print(f"警告：计算的反应进行程度不合理: x = {x_solution}")
    # 尝试解析方程: K_eq = x^2 / (0.001-x)^2
    # sqrt(K_eq) = x / (0.001-x)
    # sqrt(K_eq) * (0.001-x) = x
    # sqrt(K_eq) * 0.001 - sqrt(K_eq) * x = x
    # sqrt(K_eq) * 0.001 = x + sqrt(K_eq) * x = x * (1 + sqrt(K_eq))
    # x = sqrt(K_eq) * 0.001 / (1 + sqrt(K_eq))
    sqrt_keq = math.sqrt(K_eq)
    x_solution = sqrt_keq * 0.001 / (1 + sqrt_keq)
    print(f"解析解: x = {x_solution:.6f} M")
else:
    print(f"数值解: 反应进行程度 x = {x_solution:.6f} M")
    
# 计算平衡时各组分浓度
atp_eq = 0.001 - x_solution
glucose_eq = 0.001 - x_solution
g6p_eq = x_solution
adp_eq = x_solution

print(f"\n平衡时各组分浓度:")
print(f"[ATP] = {atp_eq:.6f} M")
print(f"[葡萄糖] = {glucose_eq:.6f} M")
print(f"[葡萄糖-6-磷酸] = {g6p_eq:.6f} M")
print(f"[ADP] = {adp_eq:.6f} M")

# 验证平衡常数
if atp_eq > 0 and glucose_eq > 0:
    calculated_keq = (g6p_eq * adp_eq) / (atp_eq * glucose_eq)
    print(f"\n验证: 计算的平衡常数 = {calculated_keq:.2e}")
    print(f"原始平衡常数 = {K_eq:.2e}")
    print(f"误差: {abs(calculated_keq - K_eq) / K_eq * 100:.4f}%")
else:
    print("无法验证，因为某些反应物浓度接近0")

print(f"\n结论：在298.15 K, pH 9 和 0.25 M 离子强度下，")
print(f"对于反应 'ATP + 葡萄糖 = 葡萄糖-6-磷酸 + ADP'，")
print(f"ATP和葡萄糖的初始浓度均为 0.001 M 时，")
print(f"最终的平衡组分浓度为：")
print(f"ATP: {atp_eq:.6f} M")
print(f"葡萄糖: {glucose_eq:.6f} M")
print(f"葡萄糖-6-磷酸: {g6p_eq:.6f} M")
print(f"ADP: {adp_eq:.6f} M")