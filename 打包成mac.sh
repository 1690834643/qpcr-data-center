#!/bin/bash
# 在 macOS 上把 qPCR 数据中心打包成 .app 并封装为 dmg。
# 用法：bash 打包成mac.sh [产物后缀]，例如 bash 打包成mac.sh arm64
set -e
SRC="$(cd "$(dirname "$0")" && pwd)"
SUFFIX="${1:-$(uname -m)}"
BUILD="$(mktemp -d /tmp/qpcr_center_build.XXXX)"
mkdir -p "$BUILD/web"
cp "$SRC"/qpcr_core.py "$SRC"/qpcr_analysis.py "$SRC"/qpcr_store.py "$SRC"/qpcr_plots.py "$SRC"/qpcr_server.py "$BUILD/"
cp "$SRC/web/index.html" "$BUILD/web/"
cd "$BUILD"

python3 -m pip install -U pyinstaller
python3 -m PyInstaller --windowed --noconfirm --name "qPCR数据中心" \
  --osx-bundle-identifier site.lyzlab.qpcr-data-center --add-data "web:web" qpcr_server.py

APP="dist/qPCR数据中心.app"
# 后台服务不需要 Dock 图标，网页全部关闭后程序自动退出
/usr/libexec/PlistBuddy -c "Add :LSUIElement bool true" "$APP/Contents/Info.plist"
codesign --force --deep -s - "$APP"

mkdir dmg && cp -R "$APP" dmg/ && ln -s /Applications dmg/Applications
hdiutil create -volname "qPCR Data Center" -srcfolder dmg -ov -format UDZO "qPCR-Data-Center-macOS-$SUFFIX.dmg"
cp "qPCR-Data-Center-macOS-$SUFFIX.dmg" "$SRC/"
echo "完成：$SRC/qPCR-Data-Center-macOS-$SUFFIX.dmg"
