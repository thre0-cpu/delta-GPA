# 生成pMg对ATP平衡浓度影响的图表
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 读取数据
df = pd.read_csv('atp_concentration_vs_pMg.csv', encoding='utf-8')

print(f'数据已读取，共{len(df)}行')
print('列名:', df.columns.tolist())
print('\n前5行数据:')
print(df.head())

# 创建图表
plt.figure(figsize=(12, 8))
plt.plot(df['pMg'], df['ATP平衡浓度 (M)'], 'b-o', markersize=4, linewidth=1.5)
plt.xlabel('pMg', fontsize=12)
plt.ylabel('ATP平衡浓度 (M)', fontsize=12)
plt.title('不同pMg下ATP的平衡浓度 (298.15 K, pH 9, 0.25 M 离子强度)', fontsize=14)
plt.grid(True, alpha=0.3)
plt.yscale('log')  # 使用对数刻度更好地显示浓度变化

# 标注关键点
min_idx = df['ATP平衡浓度 (M)'].idxmin()
max_idx = df['ATP平衡浓度 (M)'].idxmax()

plt.annotate(f'Min: ({df.loc[min_idx, "pMg"]:.1f}, {df.loc[min_idx, "ATP平衡浓度 (M)"]:.2e})',
             xy=(df.loc[min_idx, 'pMg'], df.loc[min_idx, 'ATP平衡浓度 (M)']),
             xytext=(df.loc[min_idx, 'pMg']+1.5, df.loc[min_idx, 'ATP平衡浓度 (M)']*2),
             arrowprops=dict(arrowstyle='->', color='red', lw=0.8))

plt.annotate(f'Max: ({df.loc[max_idx, "pMg"]:.1f}, {df.loc[max_idx, "ATP平衡浓度 (M)"]:.2e})',
             xy=(df.loc[max_idx, 'pMg'], df.loc[max_idx, 'ATP平衡浓度 (M)']),
             xytext=(df.loc[max_idx, 'pMg']-2, df.loc[max_idx, 'ATP平衡浓度 (M)']/2),
             arrowprops=dict(arrowstyle='->', color='red', lw=0.8))

# 绘制生理条件标注线
physiological_pmg_idx = (df['pMg'] - 3.0).abs().idxmin()
plt.axvline(x=3.0, color='green', linestyle='--', alpha=0.6, label=f'生理pMg={df.loc[physiological_pmg_idx, "pMg"]:.1f}')
plt.legend()

plt.tight_layout()
plt.savefig('atp_concentration_vs_pMg_plot.png', dpi=300, bbox_inches='tight')
print('\n图表已保存为atp_concentration_vs_pMg_plot.png')

plt.show()

# 输出关键信息
print(f'\n关键信息:')
print(f'生理条件(pMg=3)时ATP浓度: {df.loc[physiological_pmg_idx, "ATP平衡浓度 (M)"]:.8f} M')
print(f'最低ATP浓度: {df["ATP平衡浓度 (M)"].min():.8f} M (在pMg={df.loc[df["ATP平衡浓度 (M)"].idxmin(), "pMg"]}处)')
print(f'最高ATP浓度: {df["ATP平衡浓度 (M)"].max():.8f} M (在pMg={df.loc[df["ATP平衡浓度 (M)"].idxmax(), "pMg"]}处)')