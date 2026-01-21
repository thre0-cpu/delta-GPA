from equilibrator_api import ComponentContribution, Q_
import numpy as np
from scipy.optimize import root

# 初始化CC方法并设置条件
# 条件: 298.15 K, pH 9, 离子强度 0.25 M
cc = ComponentContribution()
cc.p_h = Q_(9.0)  # pH
cc.ionic_strength = Q_('0.25 M')  # 离子强度
cc.temperature = Q_('298.15 K')  # 温度
cc.p_mg = Q_(7.0)  # pMg 默认值

print("CC方法已初始化，条件设置完成：")
print(f"pH: {cc.p_h}")
print(f"离子强度: {cc.ionic_strength}")
print(f"温度: {cc.temperature}")
print(f"pMg: {cc.p_mg}")

# 定义反应 ATP + 葡萄糖 = 葡萄糖-6-磷酸 + ADP
# 使用通用名称尝试解析
try:
    reaction_str = "ATP + D-glucose = D-glucose-6-phosphate + ADP"
    print(f"定义反应: {reaction_str}")
    parsed_reaction = cc.parse_reaction_formula(reaction_str)
    print(f"反应解析成功: {parsed_reaction}")
except Exception as e:
    print(f"解析反应失败: {e}")
    # 如果通用名称不支持，尝试使用KEGG ID
    try:
        reaction_str = "C00002 + C00031 = C00636 + C00008"  # 使用KEGG ID
        print(f"尝试使用KEGG ID: {reaction_str}")
        parsed_reaction = cc.parse_reaction_formula(reaction_str)
        print(f"使用KEGG ID的反应解析成功: {parsed_reaction}")
    except Exception as e2:
        print(f"使用KEGG ID的反应也解析失败: {e2}")
        # 最后尝试更明确的格式
        reaction_str = "kegg:C00002 + kegg:C00031 = kegg:C00636 + kegg:C00008"
        print(f"尝试使用完整的KEGG ID格式: {reaction_str}")
        parsed_reaction = cc.parse_reaction_formula(reaction_str)
        print(f"使用完整KEGG ID格式的反应解析成功: {parsed_reaction}")

# 计算反应的ΔG'°
try:
    dg_prime = cc.standard_dg_prime(parsed_reaction)
    print(f"标准反应自由能变化 ΔG'° = {dg_prime}")
    
    # 计算平衡常数
    RT = cc.RT
    keq = np.exp(-dg_prime.value.m / RT.m)
    print(f"平衡常数 Keq = {keq}")
except Exception as e:
    print(f"计算反应自由能失败: {e}")

# 设置初始浓度
atp_init = Q_('0.001 M')
glucose_init = Q_('0.001 M')
g6p_init = Q_('0.0 M')
adp_init = Q_('0.0 M')

print(f"初始浓度:")
print(f"[ATP] = {atp_init}")
print(f"[葡萄糖] = {glucose_init}")
print(f"[葡萄糖-6-磷酸] = {g6p_init}")
print(f"[ADP] = {adp_init}")

# 计算平衡时的浓度
# 对于反应 aA + bB ⇌ cC + dD，平衡常数 K = ([C]^c [D]^d)/([A]^a [B]^b)
# 反应: ATP + 葡萄糖 ⇌ 葡萄糖-6-磷酸 + ADP
# 设反应进行程度为 x，则:
# [ATP] = 0.001 - x
# [葡萄糖] = 0.001 - x
# [G6P] = x
# [ADP] = x
# Keq = [G6P][ADP] / ([ATP][葡萄糖]) = x^2 / ((0.001-x)^2)

# 对于这个方程: x^2 / ((0.001-x)^2) = Keq
# 我们可以重写为: x / (0.001-x) = sqrt(Keq)
# 进一步得到: x = sqrt(Keq) * (0.001-x) = sqrt(Keq) * 0.001 - sqrt(Keq) * x
# 合并得: x + sqrt(Keq) * x = sqrt(Keq) * 0.001
# 即: x * (1 + sqrt(Keq)) = sqrt(Keq) * 0.001
# 所以: x = (sqrt(Keq) * 0.001) / (1 + sqrt(Keq))

import math

sqrt_keq = math.sqrt(keq)
x = (sqrt_keq * atp_init.m) / (1 + sqrt_keq)

print(f"反应进行程度 x = {x} M")

# 由于ATP和葡萄糖的初始浓度相同，我们可以检查x是否超过初始浓度
if x > atp_init.m:
    print("注意：反应进行程度超过了初始浓度，说明几乎完全反应")
    x = atp_init.m  # 限制反应程度不超过初始浓度

# 计算平衡浓度
atp_eq = atp_init.m - x
glucose_eq = glucose_init.m - x
g6p_eq = g6p_init.m + x
adp_eq = adp_init.m + x

print(f"平衡浓度:")
print(f"[ATP] = {atp_eq} M")
print(f"[葡萄糖] = {glucose_eq} M")
print(f"[葡萄糖-6-磷酸] = {g6p_eq} M")
print(f"[ADP] = {adp_eq} M")

# 验证平衡常数
calculated_keq = (g6p_eq * adp_eq) / (atp_eq * glucose_eq)
print(f"验证计算的平衡常数: {calculated_keq}")
print(f"从ΔG'°计算的平衡常数: {keq}")
print(f"两者是否接近: {abs(calculated_keq - keq) < 1e-6}")

# 输出最终结果
print("最终平衡组分:")
print(f"ATP: {atp_eq:.6f} M")
print(f"葡萄糖: {glucose_eq:.6f} M")
print(f"葡萄糖-6-磷酸: {g6p_eq:.6f} M")
print(f"ADP: {adp_eq:.6f} M")

# 计算反应进行程度百分比
if atp_init.m > 0:
    progress_percent = (x / atp_init.m) * 100
    print(f"\n反应进行程度: {progress_percent:.2f}%")
else:
    print(f"\n反应进行程度无法计算（初始浓度为0）")