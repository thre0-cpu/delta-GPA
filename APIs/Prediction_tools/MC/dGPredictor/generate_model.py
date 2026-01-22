"""
生成 dGPredictor 模型文件
"""
import os
import sys

# 设置路径
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, 'data')
MODEL_DIR = os.path.join(SCRIPT_DIR, 'model')

# 确保 model 目录存在
os.makedirs(MODEL_DIR, exist_ok=True)

from scipy.io import loadmat
import numpy as np
from sklearn.linear_model import BayesianRidge
from sklearn.metrics import mean_squared_error, r2_score
import joblib

print("加载训练数据...")
ac = loadmat(os.path.join(DATA_DIR, 'Test_KEGG_all_grp.mat'))

y = ac['y'].flatten()
X = ac['X_comb_all']

print(f"训练数据: X shape = {X.shape}, y shape = {y.shape}")

print("训练贝叶斯岭回归模型...")
model = BayesianRidge(tol=1e-6, fit_intercept=False, compute_score=True)
model.fit(X, y)

# 评估
y_pred = model.predict(X)
mse = mean_squared_error(y, y_pred)
r2 = r2_score(y, y_pred)

print(f"训练结果:")
print(f"  MSE = {mse:.2f}")
print(f"  R² = {r2:.4f}")

# 保存模型
model_path = os.path.join(MODEL_DIR, 'M12_model_BR.pkl')
joblib.dump(model, model_path, compress=3)
print(f"模型已保存至: {model_path}")
