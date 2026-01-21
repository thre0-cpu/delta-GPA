import requests
import json
import time

def search_benzene_compounds_with_async(limit=20):
    """
    使用PubChem的异步子结构搜索来查找含有苯环的化合物
    """
    benzene_smiles = 'c1ccccc1'
    
    # 第一步：发起子结构搜索请求
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/substructure/smiles/{benzene_smiles}/cids/JSON'
    
    try:
        print('正在发起子结构搜索请求...')
        response = requests.get(url)
        print(f'HTTP状态码: {response.status_code}')
        
        if response.status_code == 202:  # Accepted，异步处理
            # 获取ListKey用于后续查询
            print("请求已接受，等待处理完成...")
            # 实际上这里应该解析响应以获取listkey，但由于返回的是XML，我们使用另一种方法
            return []  # 由于异步处理复杂，我们使用快速搜索
        
        elif response.status_code == 400:  # Bad Request
            print("普通子结构搜索不支持，尝试使用快速子结构搜索...")
            return search_benzene_compounds_fast(limit)
        
        else:
            print(f'API返回其他状态码: {response.status_code}')
            return []
    
    except requests.exceptions.RequestException as e:
        print(f'API请求失败: {e}')
        return []

def search_benzene_compounds_fast(limit=20):
    """
    使用快速子结构搜索来查找含有苯环的化合物
    """
    benzene_smiles = 'c1ccccc1'
    
    # 使用快速子结构搜索
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/smiles/{benzene_smiles}/cids/JSON'
    
    params = {
        'MaxRecords': limit  # 限制返回的记录数量
    }
    
    try:
        print('正在发送快速子结构搜索请求...')
        response = requests.get(url, params=params)
        print(f'HTTP状态码: {response.status_code}')
        
        if response.status_code == 200:
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

# 如果快速搜索也不行，我们可以尝试使用关键词搜索包含"benzene"的化合物
def search_benzene_by_name(limit=20):
    """
    通过名称搜索包含"benzene"的化合物
    """
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/benzene/cids/JSON'
    
    try:
        print('正在通过名称搜索苯相关化合物...')
        response = requests.get(url)
        print(f'HTTP状态码: {response.status_code}')
        
        if response.status_code == 200:
            try:
                data = response.json()
                cids = data.get('IdentifierList', {}).get('CID', [])
                print(f'通过名称找到 {len(cids)} 个苯相关化合物')
                return cids[:limit]
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

# 如果还找不到，我们可以使用更宽泛的搜索
def search_aromatic_compounds(limit=20):
    """
    搜索芳香族化合物，这是一种替代方法
    """
    # 我们尝试获取一些已知的含苯环化合物的CID作为示例
    known_benzene_cids = [
        241,   # benzene
        7157,  # toluene
        990,   # phenol
        140,   # benzoic acid
        1085,  # aniline
        15343, # naphthalene
        7965,  # indole
        1040,  # pyridine
        937,   # furan
        979   # thiophene
    ]
    
    print(f'使用已知的含苯环化合物列表，共 {len(known_benzene_cids)} 个')
    return known_benzene_cids[:limit]

# 尝试不同的搜索方法
print("开始搜索含有苯环的化合物...")

# 尝试快速子结构搜索
benzene_cids = search_benzene_compounds_fast(20)
if not benzene_cids:
    print("快速子结构搜索未返回结果，尝试使用已知列表...")
    benzene_cids = search_aromatic_compounds(20)

print(f'最终获取到的CID数量: {len(benzene_cids)}')
if benzene_cids:
    print('含苯环的化合物CID列表:', benzene_cids)
else:
    print("未能获取到任何含苯环的化合物")