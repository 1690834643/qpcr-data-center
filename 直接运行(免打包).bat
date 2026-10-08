@echo off
chcp 65001 >nul
title qPCR 数据中心
cd /d "%~dp0"
REM 本机已装 Python 3 时双击运行，会自动打开浏览器。关闭此窗口即退出程序。
python qpcr_server.py
if errorlevel 1 (
  echo.
  echo 启动失败。请确认已安装 Python 3.x（安装时勾选 Add to PATH）。
  pause
)
