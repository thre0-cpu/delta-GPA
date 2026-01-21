@echo off

set "PROJECT_DIR=D:\threo333\Projects\AAA_delta_GPA\codes\qwen-web-local"
set "CONDA_BASE=D:\Anaconda"
set "CONDA_ENV=dgpa"

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
if not exist "%CONDA_BASE%\Scripts\activate.bat" (
    echo [ERROR] Conda activate.bat not found at: %CONDA_BASE%\Scripts\activate.bat
    pause
    exit /b 1
)

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

REM 添加 npm 到 PATH
set "PATH=C:\Users\threo\AppData\Roaming\npm;%PATH%"

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
