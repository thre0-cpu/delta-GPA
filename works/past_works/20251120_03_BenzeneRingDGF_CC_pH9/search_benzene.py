import requests
import json
import time

# 使用子结构搜索来查找含有苯环的化合物
def search_benzene_compounds(limit=20):
    benzene_smiles = 'c1ccccc1'
    
    # 使用substructure搜索
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/substructure/smiles/{benzene_smiles}/cids/JSON'
    
    params = {
        'MaxRecords': limit,  # 限制返回的记录数量
    }
    
    try:
        print('正在发送请求到PubChem...')
        response = requests.get(url, params=params)
        print(f'HTTP状态码: {response.status_code}')
        
        if response.status_code == 200:
            # 尝试解析JSON响应
            try:
                data = response.json()
                cids = data.get('IdentifierList', {}).get('CID', [])
                print(f'找到 {len(cids)} 个含有苯环的化合物')
                return cids[:limit]  # 返回前limit个ID
            except json.JSONDecodeError:
                print('响应不是有效的JSON格式')
                print(f'响应内容: {response.text[:500]}...')
                return []
        else:
            print(f'API返回错误状态码: {response.status_code}')
            print(f'响应内容: {response.text[:500]}...')
            return []
        
    except requests.exceptions.RequestException as e:
        print(f'API请求失败: {e}')
        return []
    except Exception as e:
        print(f'发生其他错误: {e}')
        return []

# 搜索含有苯环的化合物
print("开始搜索含有苯环的化合物...")
benzene_cids = search_benzene_compounds(20)
print(f'最终获取到的CID数量: {len(benzene_cids)}')
if benzene_cids:
    print('前20个含有苯环的化合物CID:', benzene_cids)
else:
    print("未能获取到任何含有苯环的化合物")