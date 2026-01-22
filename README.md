# delta-GPA

面向**生化热力学**的可运行 **Agent/服务** 项目。  
本仓库包含 prompts、工具接口（APIs），以及一个本地 Web 前端与服务端（位于 `qwen-web-local/`），通过 **Windows 脚本 `start.bat`** 一键构建并启动。

---

## 依赖

### 1) Conda + 指定环境
`qwen-web-local/start.bat` 会：
- 初始化 Conda
- 激活环境：`dgpa`

你需要：
- 已安装 Anaconda/Miniconda
- 已创建 conda 环境 `dgpa`

> 环境依赖以仓库内 `environment.yml` / `requirements.txt` 为准。

### 2) Node.js + npm
脚本会执行：
- `node --version`
- `npm run build`
- `npx tsx server.ts`

你需要：
- **Node.js**（建议 20+）
- npm 可用
- 项目依赖已安装（通常需要先 `npm install`）

> 脚本里还会把用户级 npm bin 目录加进 PATH：  
`C:\Users\threo\AppData\Roaming\npm`

### 3) Qwen Code（用于 Agent 能力）
Qwen Code 仓库：<https://github.com/QwenLM/qwen-code>

建议全局安装（npm）：
````bash
npm install -g @qwen-code/qwen-code@latest
````
并在首次使用时通过 `qwen` + `/auth` 完成认证（或用 OpenAI 兼容环境变量方式）。

---

## 快速开始（Windows）

### 0) 克隆仓库
````bash
git clone https://github.com/thre0-cpu/delta-GPA.git
cd delta-GPA
````

### 1) 准备 Python/Conda 环境
创建并激活 `dgpa`（示例）：
````bash
conda env create -f environment.yml
conda activate dgpa
````

### 2) 安装前端/服务端依赖（qwen-web-local）
进入目录并安装依赖：
````bash
cd qwen-web-local
npm install
````

### 3) 启动（构建前端 + 启动服务端）
直接运行：
````bat
cd qwen-web-local
start.bat
````

脚本会自动：
1. 切换到 `PROJECT_DIR`
2. 初始化 conda 并激活 `dgpa`
3. `npm run build` 构建前端
4. `npx tsx server.ts` 启动服务端

启动成功后，访问地址：**http://localhost:915**

---

## start.bat 说明

脚本会**自动检测**以下路径：

- `PROJECT_DIR`：自动获取脚本所在目录
- `CONDA_BASE`：自动检测常见 Conda 安装位置（Anaconda3、miniconda3 等）
- `CONDA_ENV`：默认为 `dgpa`（可在脚本中修改）

如果自动检测失败，脚本会提示你手动设置 `CONDA_BASE` 环境变量

---

## 目录结构

- `APIs/`：工具与接口封装（热力学计算/查询/外部服务对接）
- `envs/`：Python 环境与自定义包（如 GNN 模型）
- `prompts/`：Agent prompts 与模板
- `qwen-web-local/`：本地 Web 前端 + Node/TS 服务端（`server.ts`）+ `start.bat`
- `sl/`：辅助文档与说明
- `works/`：实验/任务产出工作区
- `environment.yml`：Conda 环境配置
- `requirements.txt`：Python 依赖列表

---

## 常见问题（FAQ）

### 1) `Conda activate.bat not found`
确认 `CONDA_BASE` 指向你的 Anaconda/Miniconda 安装目录，且存在：
`%CONDA_BASE%\Scripts\activate.bat`

### 2) `npm` / `node` 找不到
- 先安装 Node.js
- 确认 `node --version` 在命令行可用
- 如使用用户级全局安装，确保 npm bin 在 PATH（脚本已添加了示例路径）

### 3) `npm run build` 失败
通常是未安装依赖或 Node 版本不兼容：
````bash
cd qwen-web-local
npm install
npm run build
````

### 4) `npx tsx server.ts` 启动失败
- 确认 `tsx` 相关依赖已在 `package.json` 中声明并已安装
- 查看控制台报错堆栈定位缺失的环境变量/端口冲突等问题

---

## License

This project is licensed under the MIT License. See the LICENSE file for details.