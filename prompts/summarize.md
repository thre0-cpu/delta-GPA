# delta_GPA 总结任务

本文件有关delta_GPA如何总结一个项目，当用户要求**根据此文档指示进行任务总结**时，delta_GPA需要按照下面的步骤进行总结。

## 工作流与交付 (Workflow & Deliverables)

*本部分规定任务结束时的标准化动作。*

  * **任务收尾清单**：在用户要求总结整个项目时，在TODO中添加以下条目并依次执行
      * 新建一个文件夹`temp_files`，存放所有**临时脚本、Notebook**等文件；
      * 新建一个文件夹`results`，存放所有**表格、图片**等重要结果文件；
      * 新建一个文件夹`reports`；
      * 将任务的**全部工作流程**、**有效代码**和**图片**等核心文件整理到`main.ipynb`这个Jupyter Notebook中，如果没有找到则新建一个；
      * 在`reports`文件夹中创建技术报告，命名为`TECH_report.md`；
      * 在`reports`文件夹中创建代码分析报告，命名为`CODE_report.md` ；
      * 创建工作区`README.md`；
      * 确认最外层只有`main.ipynb`和`README.md`，不要在外面留下任何python脚本等文件

**注意**：整理后执行`jupyter`命令检查`main.ipynb`能否正常打开为Jupyter Notebook，特别注意中文符号和斜杠的问题。

## 项目结构

整理后的最终项目结构如下：

```
├── main.ipynb                 # 核心工作流程和可视化
├── temp_files/                # 临时脚本/Notebook文件夹
│   ├── XXX.py
│   ├── YYY.ipynb
│   └── ...
├── results/                   # 结果文件夹
│   ├── XXX.csv
│   ├── YYY.png
│   └── ...
├── reports/                   # 结果文件夹
│   ├── TECH_report.md         # 技术报告
│   └── CODE_report.md         # 代码分析报告
└── README.md                  # 说明文件
```

**注意**：最外层只有`main.ipynb`和`README.md`，不要在外面留下任何python脚本等文件

## ⚠️ 严禁幻觉执行 (CRITICAL: No Hallucinated Actions)

**这是最高优先级规则，违反此规则会导致严重后果。**

1. **说到必须做到**：当你声称"创建了文件"、"执行了代码"、"查看了文档"时，你**必须实际调用相应的工具**。绝对禁止只在文字中描述行动而不真正执行。

2. **禁止想象性完成**：以下行为是严格禁止的：
   - ❌ 说"我已经创建了 xxx.ipynb 文件"但没有调用 `create_file` 工具
   - ❌ 说"我查看了文档"但没有调用 `read_file` 工具
   - ❌ 说"代码已执行成功"但没有调用 `run_notebook_cell` 或 `run_in_terminal`
   - ❌ 列出"已完成的任务"清单但实际没有任何工具调用

3. **先做后说**：正确的流程是：
   - ✅ 先调用工具执行操作
   - ✅ 等待工具返回结果
   - ✅ 再向用户报告完成情况

4. **自检规则**：在声称完成任务之前，检查自己是否真的调用了工具。如果没有工具调用记录，就不能声称任务已完成。