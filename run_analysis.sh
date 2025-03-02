#!/bin/bash

# 檢查參數數量
if [ "$#" -ne 10 ]; then
    echo "錯誤: 參數不完整"
    echo "用法: $0 -n <normal.bam> -t <tumor.bam> -v <somatic.vcf> -r <reference.fa> -o <output.txt>"
    exit 1
fi

# 初始化變數
NORMAL_BAM=""
TUMOR_BAM=""
VCF_FILE=""
REFERENCE=""
OUTPUT_FILE=""

# 解析參數
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -n) NORMAL_BAM="$2"; shift 2 ;;
        -t) TUMOR_BAM="$2"; shift 2 ;;
        -v) VCF_FILE="$2"; shift 2 ;;
        -r) REFERENCE="$2"; shift 2 ;;
        -o) OUTPUT_FILE="$2"; shift 2 ;;
        *) echo "未知的參數: $1"; exit 1 ;;
    esac
done

# 檢查檔案是否存在
for file in "$NORMAL_BAM" "$TUMOR_BAM" "$VCF_FILE" "$REFERENCE"; do
    if [[ ! -f "$file" ]]; then
        echo "錯誤: 檔案 $file 不存在"
        exit 1
    fi
done

# 檢查主程式是否可執行
if [[ ! -x "./somatic_methylation_analyzer" ]]; then
    echo "錯誤: ./somatic_methylation_analyzer 無法執行，請檢查權限"
    exit 1
fi

# 執行程式並記錄輸出
./somatic_methylation_analyzer -o "$OUTPUT_FILE" -n "$NORMAL_BAM" -t "$TUMOR_BAM" -v "$VCF_FILE" -r "$REFERENCE" 2>&1 | tee run_analysis.log

# 檢查執行結果
if [ $? -eq 0 ]; then
    echo "分析完成，結果輸出至 $OUTPUT_FILE"
else
    echo "分析失敗，詳細錯誤請查看 run_analysis.log"
    exit 1
fi
