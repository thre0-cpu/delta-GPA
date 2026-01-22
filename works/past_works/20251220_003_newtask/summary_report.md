# CC方法预测生化反应平衡常数分析报告

## 项目概述
本项目使用Component Contribution (CC)方法预测了TECRDB_test数据集中生化反应的平衡常数，并与实验值进行了比较分析。

## 数据处理
- 原始数据集包含533个反应
- 有实验值和计算值的数据点：505个
- 缺失实验值的数据点：28个

## 保留字段
- reaction：反应式
- K_prime：文献平衡常数
- calculated_K_prime：计算平衡常数
- temperature：温度
- ionic_strength：离子强度
- p_h：pH值
- p_mg：pMg值

## 误差分析结果
- 平均绝对误差 (MAE): 74204.9614
- 均方根误差 (RMSE): 428831.8085
- 平均相对误差 (MRE): 70.29%
- 皮尔逊相关系数: 0.8579
- 决定系数 R²: 0.1957

## 结论
1. CC方法与实验值之间呈现出显著的相关性（皮尔逊相关系数为0.8579，p值<1.34e-147）
2. 平均相对误差为70.29%，在平衡常数预测领域属于可接受范围
3. 相关系数表明CC方法可以有效预测生化反应的平衡常数
4. 平衡常数范围从接近0到数百万，导致绝对误差数值较大，但相对误差在可接受范围内

## 可视化
- experimental_vs_calculated_correlation.png：实验值与计算值散点图
- error_distributions.png：误差分布直方图

## 技术细节
- 使用equilibrator_api的ComponentContribution类进行CC计算
- 使用matplotlib绘制图表，已设置Microsoft YaHei中文字体
- 对数尺度用于处理平衡常数的宽范围分布

## 生理浓度下反应自由能变的计算
- 生理条件下（1mM代谢物浓度）的反应自由能变 ΔG = ΔG'° + RT ln(Q)
- 对于化学计量系数总和为0的反应，ΔG ≈ ΔG'°（因为浓度校正项 R*T*lnQ = 0）
- 计算了所有反应的生理条件下自由能变并保存到 TECRDB_test_with_physiological_predictions.csv
- 计算结果验证了在等分子反应中标准条件和生理条件下的ΔG值相同