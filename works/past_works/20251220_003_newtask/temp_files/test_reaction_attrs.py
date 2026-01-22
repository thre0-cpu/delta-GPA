#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
测试equilibrator_api中Reaction对象的属性
"""

from equilibrator_api import ComponentContribution, Q_

def test_reaction_properties():
    print('正在初始化ComponentContribution...')
    cc = ComponentContribution()
    print('ComponentContribution初始化完成')
    
    # 解析一个简单反应来检查其属性
    rxn_str = "kegg:C00002 + kegg:C00001 = kegg:C00008 + kegg:C00009"  # ATP + H2O = ADP + Pi
    reaction = cc.parse_reaction_formula(rxn_str)
    
    print(f'反应: {rxn_str}')
    print(f'Reaction对象类型: {type(reaction)}')
    print(f'Reaction对象的所有属性和方法:')
    
    # 打印所有属性
    for attr in dir(reaction):
        if not attr.startswith('_'):  # 忽略私有属性
            try:
                value = getattr(reaction, attr)
                if not callable(value):  # 只显示非方法属性
                    print(f'  {attr}: {value}')
            except:
                print(f'  {attr}: <无法获取>')
    
    print(f'\n化学计量系数信息:')
    # 尝试找出正确的属性名
    possible_attrs = ['reac', 'prod', 'compound_stoichiometry', 'stoichiometry', 'sparse', 'sparse_representation']
    for attr in possible_attrs:
        if hasattr(reaction, attr):
            print(f'{attr}: {getattr(reaction, attr)}')
        else:
            print(f'{attr}: 不存在')

if __name__ == '__main__':
    test_reaction_properties()