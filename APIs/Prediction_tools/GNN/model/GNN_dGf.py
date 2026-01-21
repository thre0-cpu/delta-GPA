import sys
import json
import os

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import numpy as np

import torch
import torch.nn as nn

from .dGbyG import mol_to_graph_data, MP_network, threo, best_model_params

class networks(nn.Module):
    def __init__(self, dir) -> None:
        super().__init__()
        self.nets = nn.ModuleList([])
        for file in os.listdir(dir):
            path = os.path.join(dir, file)
            net = MP_network(atom_dim=139, bond_dim=23, emb_dim=300, num_layer=2)
            net.load_state_dict(torch.load(path, map_location=torch.device('cpu')))
            self.nets.append(net)
        self.num = len(self.nets)
    
    def forward(self, data):
        outputs = torch.zeros(size=(self.num, 1)).to(data.x.device)  # shape=[number of net, 1]

        for i, net in enumerate(self.nets):
            outputs[i] = net(data)  # net.shape = [1,1] or [atom number, 1]

        return outputs  # .squeeze()


# 使用方式
network = networks(best_model_params)
network.eval()
device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
network.to(device)

def predict_standard_dGf_prime(
    inchi: str, 
) -> tuple:
    """
    预测标准吉布斯自由能
    
    Args:
        inchi: 分子的 InChI 字符串
    
    Returns:
        (mean, std)
    """
    try:
        # 1. 从 InChI 读取分子
        mol = Chem.MolFromInchi(inchi, removeHs=False, sanitize=True)
        if mol is None:
            raise ValueError(f"无法解析 InChI: {inchi}")
        
        # 2. 第一次 normalize（Normalize + Uncharger）
        mol = rdMolStandardize.Normalize(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        
        # 3. 转换成 SMILES 再读回来（关键步骤！）
        smiles = Chem.MolToSmiles(mol)
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        
        # 4. 第二次 normalize
        mol = rdMolStandardize.Normalize(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        
        # 5. 添加氢原子
        mol = Chem.AddHs(mol)
        
        # 6. 转换为图数据
        data = mol_to_graph_data(mol).to(device)
        
        # 7. 预测
        with torch.no_grad():
            predictions = network(data).cpu().numpy()
        
        mean_pred = np.mean(predictions)
        std_pred = np.std(predictions)
        return mean_pred, std_pred
        
    
    except Exception as e:
        print(f"预测失败 - InChI: {inchi}, 错误: {e}")
        return np.nan, np.nan
