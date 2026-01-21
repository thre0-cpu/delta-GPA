# Pub Chem 数据库使用指南

如果本文件夹被引用，说明用户需要使用PubChem数据库中的数据。
PubChem提供REST API接口，可以用Python的request相关命令实现调用，**无需安装其他包**。
具体的API文档和调用方式参考文件夹内另外两个文档，**一定要阅读**
不要使用**pubchempy**

## 📗 Document 1: PUG REST Tutorial（教程文档）
@ 01_PUG REST Tutorial.md
**用途**: 日常查询的首选参考
**特点**: 
- 示例丰富，易于理解
- 按使用场景组织
- 包含最新功能详解

**何时使用**:
- ✅ 需要了解如何完成某个任务
- ✅ 寻找示例URL
- ✅ 学习新功能（Gene/Protein/Pathway）
- ✅ 理解概念和工作流程

## 📕 Document 2: PUG REST Specification（规范文档）
@ 02_PUG REST Specification.md
**用途**: 技术细节和完整参考
**特点**:
- 参数定义精确
- 完整的选项列表
- 系统化的技术规范

**何时使用**:
- ✅ 需要完整的参数列表
- ✅ 查询HTTP状态码含义
- ✅ 确认语法细节
- ✅ 查找不常用的选项

---

## 🎯 查询策略

### 策略1: 快速任务（80%的情况）
```
步骤：
1. 在Tutorial中搜索相关示例
2. 根据示例构建URL
3. 执行并返回结果
```

### 策略2: 参数确认（15%的情况）
```
步骤：
1. 在Tutorial中找到基本方法
2. 在Specification中确认参数细节
3. 构建精确的URL
4. 执行并返回结果
```

### 策略3: 复杂场景（5%的情况）
```
步骤：
1. 在Tutorial中理解整体流程
2. 在Specification中查找所有相关选项
3. 设计多步骤解决方案
4. 逐步执行并处理中间结果
```

---

## 🔍 文档搜索技巧

### 按关键词搜索
```
需要查询...              搜索关键词
─────────────────────    ─────────────────
化合物属性               "property"
结构图                   "PNG" 或 "image"
同义词                   "synonyms"
交叉引用                 "xref"
异步操作                 "ListKey" 或 "asynchronous"
批量查询                 "POST" 或 "list"
错误处理                 "status code" 或 "error"
```

### 按URL模式搜索
```
如果用户提到...          搜索模式
─────────────────────    ─────────────────
"通过名称"               "/name/"
"通过SMILES"             "/smiles/"
"分子量"                 "MolecularWeight"
"子结构"                 "substructure"
"相似"                   "similarity"
"基因"                   "/gene/"
"测定"                   "/assay/"
```

---