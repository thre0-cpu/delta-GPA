#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
处理TECRDB_test_with_predictions.csv文件，只保留指定的列
"""

import pandas as pd

def main():
    # 读取预测结果文件
    df = pd.read_csv('TECRDB_test_with_predictions_full.csv')
    
    print(f'原始数据集包含 {len(df)} 行')
    print('原始列包括：', list(df.columns))
    
    # 保留指定的列
    selected_columns = ['reaction', 'K_prime', 'calculated_K_prime', 'temperature', 'ionic_strength', 'p_h', 'p_mg']
    
    # 检查所有需要的列是否存在
    missing_columns = [col for col in selected_columns if col not in df.columns]
    if missing_columns:
        print(f'警告: 以下列不存在于数据中: {missing_columns}')
        available_columns = [col for col in selected_columns if col in df.columns]
        df_filtered = df[available_columns]
    else:
        df_filtered = df[selected_columns]
    
    # 保存处理后的数据
    df_filtered.to_csv('TECRDB_test_filtered.csv', index=False)
    print(f'处理完成，保留了 {len(df_filtered)} 行数据')
    print('保留的列：', list(df_filtered.columns))
    print('前5行数据预览：')
    print(df_filtered.head())

if __name__ == '__main__':
    main()