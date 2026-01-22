#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
检查Compound对象的属性
"""

from equilibrator_api import ComponentContribution

def check_compound_attributes():
    print('正在初始化ComponentContribution...')
    cc = ComponentContribution()
    print('ComponentContribution初始化完成')
    
    # 解析一个简单反应来检查Compound对象的属性
    rxn_str = "kegg:C00002 + kegg:C00001 = kegg:C00008 + kegg:C00009"  # ATP + H2O = ADP + Pi
    reaction = cc.parse_reaction_formula(rxn_str)
    
    print(f'反应: {rxn_str}')
    print('反应中的Compound对象属性:')
    
    for compound, coeff in reaction.sparse.items():
        print(f'\\nCompound对象类型: {type(compound)}')
        print(f'化学计量系数: {coeff}')
        
        # 检查可用属性
        for attr in dir(compound):
            if not attr.startswith('_'):
                try:
                    value = getattr(compound, attr)
                    if not callable(value):
                        print(f'  {attr}: {value}')
                except:
                    print(f'  {attr}: <无法获取>')

if __name__ == '__main__':
    check_compound_attributes()