@echo off
chcp 65001 >nul
title qPCR 数据中心 - 打包 exe
echo ============================================
echo   qPCR 数据中心  打包为单个 exe
echo ============================================
echo.
REM 在英文临时目录打包，避免中文路径让 PyInstaller 出错
set BUILD=C:\qpcr_center_build
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\web"
copy /y "%~dp0qpcr_core.py" "%BUILD%\" >nul
copy /y "%~dp0qpcr_analysis.py" "%BUILD%\" >nul
copy /y "%~dp0qpcr_store.py" "%BUILD%\" >nul
copy /y "%~dp0qpcr_plots.py" "%BUILD%\" >nul
copy /y "%~dp0qpcr_server.py" "%BUILD%\" >nul
copy /y "%~dp0web\index.html" "%BUILD%\web\" >nul
cd /d "%BUILD%"

echo [1/3] 安装打包工具 PyInstaller（清华镜像）...
python -m pip install -U pyinstaller -i https://pypi.tuna.tsinghua.edu.cn/simple
if errorlevel 1 (
  echo 安装失败：请确认已安装 Python 并勾选了 Add to PATH。也可换镜像 https://mirrors.aliyun.com/pypi/simple
  if not "%1"=="nopause" pause
  exit /b 1
)

echo.
echo [2/3] 正在打包（约 1-2 分钟）...
python -m PyInstaller --onefile --noconsole --name "qPCR数据中心" --add-data "web;web" qpcr_server.py
if errorlevel 1 (
  echo 打包失败。
  if not "%1"=="nopause" pause
  exit /b 1
)

echo.
echo [3/3] 复制 exe 回原文件夹...
copy /y "%BUILD%\dist\qPCR数据中心.exe" "%~dp0" >nul
echo.
echo 完成！双击 qPCR数据中心.exe 即可使用（可单独发给同学，无需装 Python）
if not "%1"=="nopause" pause
