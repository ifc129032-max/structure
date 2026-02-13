# 混合型連接梁分析專案

此專案提供一個以 **Tkinter** 建立的桌面工具，用於「軟化壓拉桿模型」下的混合型連接梁設計計算，包含：

- 幾何尺寸、材料強度、需求參數輸入
- 詳細計算書文字輸出（內建於 GUI）
- 使用 ReportLab 匯出 PDF（支援圖文並排區塊）

## 專案結構

- `app.py`：主程式（GUI、計算核心、PDF 匯出）
- `requirements.txt`：相依套件

> 若要在 PDF 中嵌入示意圖，請將 `1.png`、`2.png` 放在執行目錄。

## 安裝

```bash
python -m venv .venv
source .venv/bin/activate  # Windows 請改用 .venv\Scripts\activate
pip install -r requirements.txt
```

## 執行

```bash
python app.py
```

## 注意事項

- 若未安裝 `reportlab`，程式仍可執行計算，但 PDF 匯出功能會顯示警告。
- 若環境沒有 Windows 字型（如 `msjh.ttc`），PDF 會退回預設字型。
