# 修复中文字体的ATP浓度vs pMg图表
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

# 解决中文字体问题
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

# 如果上述字体不可用，尝试从系统获取中文字体
try:
    # 在Windows系统中，常见的中文字体
    windows_font_path = [
        'C:/Windows/Fonts/simhei.ttf',  # 黑体
        'C:/Windows/Fonts/simsun.ttc',  # 宋体
        'C:/Windows/Fonts/msyh.ttc',   # 微软雅黑
        'C:/Windows/Fonts/msyhbd.ttc'  # 微软雅黑粗体
    ]
    
    # 尝试加载系统字体
    for font_path in windows_font_path:
        try:
            font_manager.fontManager.addfont(font_path)
            print(f"成功加载字体: {font_path}")
            # 使用新加载的字体
            plt.rcParams['font.sans-serif'] = [font_manager.FontProperties(fname=font_path).get_name()]
            break
        except:
            continue
except:
    print("无法加载特定中文字体，使用默认设置")

# 读取数据
df = pd.read_csv('atp_concentration_vs_pMg.csv', encoding='utf-8')

print(f'数据已读取，共{len(df)}行')
print('列名:', df.columns.tolist())

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

plt.annotate(f'最小值: ({df.loc[min_idx, "pMg"]:.1f}, {df.loc[min_idx, "ATP平衡浓度 (M)"]:.2e})',
             xy=(df.loc[min_idx, 'pMg'], df.loc[min_idx, 'ATP平衡浓度 (M)']),
             xytext=(df.loc[min_idx, 'pMg']+1.5, df.loc[min_idx, 'ATP平衡浓度 (M)']*2),
             arrowprops=dict(arrowstyle='->', color='red', lw=0.8))

plt.annotate(f'最大值: ({df.loc[max_idx, "pMg"]:.1f}, {df.loc[max_idx, "ATP平衡浓度 (M)"]:.2e})',
             xy=(df.loc[max_idx, 'pMg'], df.loc[max_idx, 'ATP平衡浓度 (M)']),
             xytext=(df.loc[max_idx, 'pMg']-2, df.loc[max_idx, 'ATP平衡浓度 (M)']/2),
             arrowprops=dict(arrowstyle='->', color='red', lw=0.8))

# 绘制生理条件标注线
physiological_pmg_idx = (df['pMg'] - 3.0).abs().idxmin()
plt.axvline(x=3.0, color='green', linestyle='--', alpha=0.6, label=f'生理pMg={df.loc[physiological_pmg_idx, "pMg"]:.1f}')
plt.legend()

plt.tight_layout()
plt.savefig('atp_concentration_vs_pMg_plot_chinese_font.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
print('\n图表已保存为atp_concentration_vs_pMg_plot_chinese_font.png')

plt.show()

# 输出关键信息
print(f'\n关键信息:')
print(f'生理条件(pMg=3)时ATP浓度: {df.loc[physiological_pmg_idx, "ATP平衡浓度 (M)"]:.8f} M')
print(f'最低ATP浓度: {df["ATP平衡浓度 (M)"].min():.8f} M (在pMg={df.loc[df["ATP平衡浓度 (M)"].idxmin(), "pMg"]}处)')
print(f'最高ATP浓度: {df["ATP平衡浓度 (M)"].max():.8f} M (在pMg={df.loc[df["ATP平衡浓度 (M)"].idxmax(), "pMg"]}处)')