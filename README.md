# Somatic Mutation 與 Methylation 分析工具

**版本：1.0.0**
開發者：liaoyoyo

## 簡介

本工具用於解析 Somatic VCF 檔案，並結合 Tumor BAM 檔案進行 Somatic Mutation 與 DNA 甲基化（Methylation）分析。利用 HTSlib 讀取與解析 BAM/VCF 檔案，並透過 OpenMP 平行運算以提升大規模數據分析效能。最終結果將分別輸出至文字檔，包含 mutation 統計與 methylation 分析資料。  
*注意：目前 Normal BAM 分析功能尚未實作，可依需求擴充。*

---

## 主要功能

- **VCF 檔案解析**  
  讀取並解析 Somatic VCF 檔案，提取每筆 mutation 之染色體、位置、參考與突變 allele 等資訊。

- **BAM 檔案分析**  
  解析 Tumor BAM 檔案，根據 CIGAR 與 MD 標籤映射 read 至參考座標，並利用 MM/ML 標籤取得 DNA 甲基化機率，進行 mutation 與 methylation 數據統計。

- **平行運算**  
  使用 OpenMP 平行化每個 somatic site 的分析，加速資料處理，同時確保多緒安全（每緒獨立開啟 BAM 檔案並於關鍵區段寫入全域結果）。

- **結果輸出**  
  產生兩份輸出檔：
  - **Somatic_analy.txt**：包含 mutation 統計資料（如染色體、位置、讀數與平均甲基化值）。
  - **methyl_analy.txt**：依照指定排序輸出 methylation 分析資料（包含相關 somatic 位點、突變類型與甲基化分數）。

- **計時工具**  
  利用 RAII 設計模式管理 Timer，以計算各流程（VCF 解析、BAM 分析與結果輸出）的執行時間，方便性能調校與除錯。

---

## 需求

- **編譯器**：支援 C++11 或更新版本  
- **HTSlib**：用於 BAM 與 VCF 檔案解析  
- **GNU getopt**：用於命令列參數解析  
- **OpenMP**：選用以進行平行運算（若編譯器支援）

---

## 編譯與安裝

請先確認系統中已安裝 HTSlib，並將其 include 與 lib 路徑設定正確。可使用下列命令進行編譯：

```bash
g++ -std=c++11 -fopenmp -I/path/to/htslib/include -L/path/to/htslib/lib main.cpp ArgParser.cpp Analysis.cpp VCFHandler.cpp Utility.cpp OutputHandler.cpp -o somatic_analysis -lhts
```

若有需要，可根據實際環境調整編譯選項與路徑設定。

---

## 使用說明

執行程式時，請依下列參數格式呼叫：

```bash
./somatic_analysis [options]
```

### 主要參數

| 參數                    | 說明                                                      |
|-------------------------|----------------------------------------------------------|
| `-t, --tumor <file>`    | 指定 Tumor BAM 檔案 (必填)                                |
| `-v, --vcf <file>`      | 指定 Somatic VCF 檔案 (必填)                              |
| `-n, --normal <file>`   | 指定 Normal BAM 檔案 (可選，目前尚未實作)                   |
| `-r, --ref <file>`      | 參考基因組檔案 (可選)                                      |
| `-o, --output <folder>` | 指定輸出資料夾 (預設：`./`)                                |
| `-w, --window <num>`    | 指定分析範圍 (預設：2000 bp)                               |
| `-j, --threads <num>`   | 設定執行緒數 (預設使用系統最大可用數)                       |
| `-h, --help`            | 顯示使用說明                                              |

---

## 輸出結果說明

- **Somatic_analy.txt**  
  每筆資料包含：  
  - 染色體 (chr)  
  - 位置 (POS)  
  - 參考 allele (ref)  
  - 突變 allele (alt)  
  - 參考讀數 (ref_count)  
  - 突變讀數 (alt_count)  
  - 參考平均甲基化值 (ref_methyl)  
  - 突變平均甲基化值 (alt_methyl)

- **methyl_analy.txt**  
  每筆資料包含：  
  - 甲基化所在染色體 (Methyl_Chr)  
  - 甲基化位點 (Methyl_POS)  
  - 相關 somatic mutation 位置 (Somatic_POS)  
  - 該筆 read 於 mutation 位點的 allele (Somatic_Allele)  
  - 平均甲基化分數 (Methylation_Score)

---

## 範例

假設有以下檔案：
 - Tumor BAM：`tumor.bam`
 -  Somatic VCF：`somatic.vcf`

執行命令範例如下：

```bash
./somatic_analysis -t tumor.bam -v somatic.vcf -o ./result_folder -w 2000 -j 4
```

程式將解析 VCF 與 BAM 檔案，執行平行分析後，於 `./result_folder` 輸出 `Somatic_analy.txt` 與 `methyl_analy.txt`，同時在終端機顯示各階段耗時資訊。

---

## 檔案結構

- **main.cpp**  
  程式進入點，依序呼叫參數解析、VCF 解析、BAM 分析與結果輸出。

- **ArgParser.cpp / ArgParser.hpp**  
  負責命令列參數的解析與說明顯示。

- **VCFHandler.cpp / VCFHandler.hpp**  
  負責解析 Somatic VCF 檔案，提取 mutation 相關資訊。

- **Analysis.cpp / Analysis.hpp**  
  執行主要分析流程，包括 read 映射、甲基化記錄解析與統計，同時使用 OpenMP 平行運算。

- **OutputHandler.cpp / OutputHandler.hpp**  
  負責將分析結果分別輸出至文字檔，並進行排序與格式化。

- **Utility.cpp / Utility.hpp**  
  提供計時工具 (Timer)，用以記錄各流程的執行時間。

---

## 注意事項

- **Normal BAM 分析尚未實作**：目前僅支援 Tumor BAM 的分析，如需使用 Normal BAM 分析，請自行擴充相關模組。
- **平行化設定**：程式內利用 OpenMP 實現平行計算，請確認編譯器支援 OpenMP。
- **資源管理**：程式採用 RAII 模式管理 Timer 等資源，確保各流程計時精確。

---

## 聯絡與貢獻

若您有任何建議或發現程式問題，歡迎提交 Issue 或 Pull Request，我們會持續優化此工具。

---

## 授權

本工具採用 **GPL v3 授權條款**
