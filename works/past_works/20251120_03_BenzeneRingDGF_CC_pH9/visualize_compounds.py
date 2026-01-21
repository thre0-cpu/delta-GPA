import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw, rdMolDescriptors
import numpy as np
from matplotlib.patches import Rectangle
import os

# 读取之前计算得到的CSV文件
df = pd.read_csv('benzene_ring_compounds_dgf_pH9.csv')

# 过滤出成功计算出自由能的化合物
successful_df = df[df['Dgf_prime'].notna()].copy()

if successful_df.empty:
    print("没有成功计算的化合物数据")
    exit()

print(f"成功计算的化合物数量: {len(successful_df)}")
print("前5行数据:")
print(successful_df[['Name', 'Dgf_prime', 'InChI']].head())

# 为每个化合物创建分子对象
successful_df = successful_df.reset_index(drop=True)
molecules = []
valid_indices = []

for idx, row in successful_df.iterrows():
    inchi = row['InChI']
    if pd.notna(inchi) and inchi != '':
        mol = Chem.MolFromInchi(inchi)
        if mol is not None:
            molecules.append(mol)
            valid_indices.append(idx)
        else:
            print(f"无法从InChI创建分子对象: {row['Name']} (CID: {row['CID']})")
    else:
        print(f"InChI为空: {row['Name']} (CID: {row['CID']})")

# 只保留能成功创建分子对象的行
valid_df = successful_df.loc[valid_indices].copy()
valid_df = valid_df.reset_index(drop=True)

print(f"成功创建分子对象的数量: {len(valid_df)}")

# 如果没有有效的分子对象，输出错误信息
if len(valid_df) == 0:
    print("没有有效的分子结构可以绘制")
    exit()

# 创建分子图像
mol_images = []
for i, mol in enumerate(molecules[:len(valid_df)]):
    try:
        img = Draw.MolToImage(mol, size=(200, 200))
        mol_images.append(img)
    except Exception as e:
        print(f"绘制分子图像时出错: {valid_df.iloc[i]['Name']}, 错误: {e}")
        # 创建一个空图像作为占位符
        from PIL import Image
        img = Image.new('RGB', (200, 200), color='white')
        mol_images.append(img)

# 创建图表
n_compounds = len(valid_df)
ncols = 4
nrows = int(np.ceil(n_compounds / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(16, 4*nrows))
if nrows == 1:
    axes = axes.reshape(1, -1)

# 为每个化合物绘制子图
for i in range(nrows):
    for j in range(ncols):
        idx = i * ncols + j
        if idx < len(valid_df):
            ax = axes[i, j]
            
            # 显示分子结构图
            ax.imshow(mol_images[idx])
            ax.axis('off')
            
            # 添加化合物名称和自由能值
            name = valid_df.iloc[idx]['Name']
            dgf = valid_df.iloc[idx]['Dgf_prime']
            cid = valid_df.iloc[idx]['CID']
            uncertainty = valid_df.iloc[idx]['Uncertainty']
            
            title = f"{name}\nCID: {cid}\nDgf'°: {dgf:.2f} ± {uncertainty:.2f} kJ/mol"
            ax.set_title(title, fontsize=10, pad=5)
        else:
            # 隐藏空的子图
            axes[i, j].axis('off')

plt.tight_layout()
plt.subplots_adjust(hspace=0.5, wspace=0.3)

# 保存图表
output_path = 'benzene_compounds_dgf_visualization.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"图表已保存到: {output_path}")

# 显示图表
plt.show()