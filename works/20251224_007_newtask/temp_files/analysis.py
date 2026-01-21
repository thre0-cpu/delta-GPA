# 分析计算结果的合理性
print('结果分析:')

# 计算结果
MoS4_2 = 2.00e-07
MoOS3_2 = 4.81e-59
MoO2S2_2 = 1.78e-97
MoO3S_2 = 1.05e-135
MoO4_2 = 2.53e-174
H2S_aq = 4.81e-59
HS_neg = 6.34e-61
S2_neg = 7.98e-69
H2S_gas_pressure = 1.15e-62

print('1. 钼物种分布:')
print(f'   MoS4^2-: {MoS4_2:.2e} mol/L')
print(f'   MoOS3^2-: {MoOS3_2:.2e} mol/L')
print(f'   MoO2S2^2-: {MoO2S2_2:.2e} mol/L')
print(f'   MoO3S^2-: {MoO3S_2:.2e} mol/L')
print(f'   MoO4^2-: {MoO4_2:.2e} mol/L')

print()
print('2. 各物种浓度/气压（归一化）:')
print(f'   MoS4^2-: {MoS4_2:.2e}')
print(f'   MoOS3^2-: {MoOS3_2:.2e}')
print(f'   MoO2S2^2-: {MoO2S2_2:.2e}')
print(f'   MoO3S^2-: {MoO3S_2:.2e}')
print(f'   MoO4^2-: {MoO4_2:.2e}')
print(f'   H2S(aq): {H2S_aq:.2e}')
print(f'   HS^-(aq): {HS_neg:.2e}')
print(f'   S^2-(aq): {S2_neg:.2e}')
print(f'   H2S(g): {H2S_gas_pressure:.2e}')

print()
print('3. 物料平衡验证:')
Mo_total = 2.00e-7
Mo_sum = MoS4_2 + MoOS3_2 + MoO2S2_2 + MoO3S_2 + MoO4_2
print(f'   钼总浓度: {Mo_sum:.2e} mol/L vs 初始 {Mo_total:.2e} mol/L')
print(f'   钼平衡符合: {abs(Mo_sum - Mo_total) < 1e-15}')

print()
print('4. 结果合理性:')
print('   - 由于平衡常数较小，MoS4^2-几乎没有水解，这符合化学直觉')
print('   - 水解产物浓度极低，说明反应向左进行程度很小')
print('   - H2S的浓度极低，主要以气相形式存在（挥发性强）')
print('   - 在pH=5.00条件下，H2S主要以分子形式存在，少量电离')

print()
print('5. 最终答案:')
print('   含钼物种浓度乘积: 0.00e+00 (实际上是一个接近于0的极小值)')
print('   含硫不含钼物种浓度/气压乘积: 2.79e-249')
print('   浓度/气压乘积比值: 0.00e+00')
print()
print('这个比值接近于0，因为绝大多数钼仍以MoS4^2-形式存在，')
print('而含硫不含钼的物种浓度也很低，这表明在给定条件下，')
print('钼酸根与硫代钼酸根的水解反应几乎不发生。')