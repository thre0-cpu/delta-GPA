@echo off
chcp 65001 >nul 2>&1

REM ============================================
REM 自动推导路径，无需手动配置
REM ============================================

REM 获取脚本所在目录（即 qwen-web-local）
set "PROJECT_DIR=%~dp0"
REM 去掉末尾的反斜杠
if "%PROJECT_DIR:~-1%"=="\" set "PROJECT_DIR=%PROJECT_DIR:~0,-1%"

REM Conda 环境名（可根据需要修改）
set "CONDA_ENV=dgpa"

REM 自动检测 Conda 安装路径
set "CONDA_BASE="
if exist "%USERPROFILE%\Anaconda3\Scripts\activate.bat" set "CONDA_BASE=%USERPROFILE%\Anaconda3"
if exist "%USERPROFILE%\miniconda3\Scripts\activate.bat" set "CONDA_BASE=%USERPROFILE%\miniconda3"
if exist "%LOCALAPPDATA%\Anaconda3\Scripts\activate.bat" set "CONDA_BASE=%LOCALAPPDATA%\Anaconda3"
if exist "%LOCALAPPDATA%\miniconda3\Scripts\activate.bat" set "CONDA_BASE=%LOCALAPPDATA%\miniconda3"
if exist "D:\Anaconda\Scripts\activate.bat" set "CONDA_BASE=D:\Anaconda"
if exist "C:\Anaconda3\Scripts\activate.bat" set "CONDA_BASE=C:\Anaconda3"
if exist "C:\ProgramData\Anaconda3\Scripts\activate.bat" set "CONDA_BASE=C:\ProgramData\Anaconda3"

REM 自动添加用户级 npm 到 PATH
if exist "%APPDATA%\npm" set "PATH=%APPDATA%\npm;%PATH%"

echo ============================================
echo Starting Qwen Web Local Server for Delta GPA
echo ============================================
echo.

cd /d "%PROJECT_DIR%"
if errorlevel 1 (
    echo [ERROR] Failed to change directory
    pause
    exit /b 1
)

echo [INFO] Current directory: %cd%
echo.

REM 检查 conda 是否存在
if "%CONDA_BASE%"=="" (
    echo [ERROR] 未能自动检测到 Conda 安装路径
    echo [INFO] 请手动设置 CONDA_BASE 环境变量，或在脚本中添加你的 Conda 路径
    pause
    exit /b 1
)
if not exist "%CONDA_BASE%\Scripts\activate.bat" (
    echo [ERROR] Conda activate.bat not found at: %CONDA_BASE%\Scripts\activate.bat
    pause
    exit /b 1
)

echo [INFO] Detected Conda at: %CONDA_BASE%
echo [INFO] Initializing conda...
call "%CONDA_BASE%\Scripts\activate.bat" "%CONDA_BASE%"
if errorlevel 1 (
    echo [ERROR] Failed to initialize conda
    pause
    exit /b 1
)

echo [INFO] Activating environment: %CONDA_ENV%
call conda activate %CONDA_ENV%
if errorlevel 1 (
    echo [ERROR] Failed to activate conda environment: %CONDA_ENV%
    pause
    exit /b 1
)

echo [INFO] Environment activated successfully
echo [INFO] Node version:
node --version
echo.

echo [INFO] Building frontend...
call npm run build
if errorlevel 1 (
    echo [ERROR] Failed to build frontend
    echo Showing build error details...
    pause
    exit /b 1
)
echo [INFO] Frontend build completed
echo.

echo [INFO] Starting server...
echo.

REM 运行服务器，捕获所有输出
npx tsx server.ts
if errorlevel 1 (
    echo.
    echo [ERROR] Server crashed with error code: %errorlevel%
)

echo.
echo [INFO] Server stopped
pause
