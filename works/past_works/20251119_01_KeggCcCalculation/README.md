# 项目：KEGG前100个化合物生成自由能计算

本项目旨在批量计算KEGG数据库中前100个化合物在pH=9条件下的标准生成吉布斯自由能(ΔG'f)。

## 文件结构

- `calculate_dg.py`
  - **描述**: 核心计算脚本。此脚本负责初始化`equilibrator-api`，遍历KEGG ID，执行计算，并生成最终的CSV结果文件。
  - **运行**: `python calculate_dg.py`

- `kegg_cc_ph9_results.csv`
  - **描述**: 最终的计算结果数据文件。包含了每个化合物的KEGG ID、计算出的ΔG'f (kJ/mol)、不确定度以及任何可能发生的错误信息。

- `technical_report.md`
  - **描述**: 技术报告。详细说明了项目的目标、方法、执行过程、最终结果和结论。

- `code_analysis_report.md`
  - **描述**: 代码分析报告。记录了在开发过程中遇到的技术挑战（如`.ipynb`生成失败、API返回值不一致等）、解决方案和对底层API的改进建议。

- `kegg_cc_calculation.ipynb`
  - **描述**: 最初尝试使用的Jupyter Notebook。**此文件已废弃**，由于执行时出现JSON格式错误，后被`calculate_dg.py`脚本取代。保留此文件用于追溯开发过程。

## 如何重现结果

1.  确保已正确配置 `dgpa` conda 环境，且 `equilibrator-api` 等依赖已安装。
2.  在项目根目录 (`codes/`) 下运行以下命令：
    ```shell
    python works/20251119_01_KeggCcCalculation/calculate_dg.py
    ```
3.  脚本执行完毕后，最新的结果将会覆盖`kegg_cc_ph9_results.csv`文件。
