import pandas as pd
from equilibrator_api import ComponentContribution, Q_
from tqdm import tqdm
import numpy as np
import warnings
import sys

# 将默认编码设置为UTF-8
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

warnings.filterwarnings('ignore')

print("正在初始化 ComponentContribution，请稍候...")
try:
    cc = ComponentContribution()
    print("初始化完成。")
except Exception as e:
    print(f"初始化失败: {e}")
    sys.exit(1)

# 设置用户指定的条件
cc.p_h = Q_(9.0)
cc.p_mg = Q_(3.0)  # 默认值
cc.ionic_strength = Q_('0.25M') # 默认值
cc.temperature = Q_('298.15K') # 默认值

print(f"计算条件已设置为: pH={cc.p_h.m}, pMg={cc.p_mg.m}, I={cc.ionic_strength.m}, T={cc.temperature.m}")

kegg_ids = [f"C{str(i).zfill(5)}" for i in range(1, 101)]
results = []

print(f"开始为 {len(kegg_ids)} 个化合物计算生成自由能...")

for kegg_id in tqdm(kegg_ids):
    try:
        compound = cc.get_compound(f"kegg:{kegg_id}")
        
        if compound is None:
            results.append({
                "KEGG_ID": kegg_id,
                "dGf_prime (kJ/mol)": np.nan,
                "Uncertainty (kJ/mol)": np.nan,
                "Error": "Compound not found in Equilibrator database"
            })
            continue
        
        standard_dgf_prime, sigma_fin, sigma_inf = cc.standard_dg_formation(compound)
        
        # 修正：检查返回值类型
        if hasattr(standard_dgf_prime, 'value'):
            dg_value = standard_dgf_prime.value.m_as('kJ/mol')
        else:
            dg_value = float(standard_dgf_prime) # 直接使用数值
            
        uncertainty = np.linalg.norm(sigma_fin) if sigma_fin is not None else 0.0
        
        results.append({
            "KEGG_ID": kegg_id,
            "dGf_prime (kJ/mol)": dg_value,
            "Uncertainty (kJ/mol)": uncertainty,
            "Error": ""
        })

    except Exception as e:
        results.append({
            "KEGG_ID": kegg_id,
            "dGf_prime (kJ/mol)": np.nan,
            "Uncertainty (kJ/mol)": np.nan,
            "Error": str(e)
        })

print("计算完成。")

results_df = pd.DataFrame(results)
output_path = 'works/20251119_01_KeggCcCalculation/kegg_cc_ph9_results.csv'
results_df.to_csv(output_path, index=False, encoding='utf-8-sig')

print(f"结果已保存到: {output_path}")
print("\n前5行结果预览:")
print(results_df.head())