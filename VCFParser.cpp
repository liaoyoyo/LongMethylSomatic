#include "VCFParser.hpp"
#include <iostream>
#include <htslib/vcf.h> // 確保包含 htslib 的 VCF 相關頭文件

// VCFParser 構造函數
VCFParser::VCFParser(const std::string &filename) : filename(filename) {
    if (filename.empty()) {
        std::cerr << "檔案名稱無效\n";
        throw std::invalid_argument("檔案名稱無效");
    }
}

// 使用 `htslib` 讀取 VCF
bool VCFParser::loadVariants() {
    htsFile *vcfFile = hts_open(filename.c_str(), "r");  
    if (!vcfFile) {
        std::cerr << "無法開啟 VCF 檔案: " << filename << "\n";
        return false; // 如果無法開啟文件，返回 false
    }

    bcf_hdr_t *hdr = bcf_hdr_read(vcfFile); 
    if (!hdr) {
        std::cerr << "無法讀取 VCF header: " << filename << "\n";
        hts_close(vcfFile); // 關閉文件
        return false; // 如果無法讀取標頭，返回 false
    }

    bcf1_t *record = bcf_init(); 
    if (!record) {
        std::cerr << "無法初始化 VCF 讀取\n";
        bcf_hdr_destroy(hdr); // 釋放標頭
        hts_close(vcfFile); // 關閉文件
        return false; // 如果無法初始化，返回 false
    }

    // 讀取 VCF 變異資料
    while (bcf_read(vcfFile, hdr, record) >= 0) {
        bcf_unpack(record, BCF_UN_ALL);  // 解碼變異資料

        // 獲取染色體名稱
        const char *chrom = bcf_hdr_id2name(hdr, record->rid);
        if (!chrom) {
            std::cerr << "無法獲取染色體名稱\n";
            continue; // 如果染色體名稱無效，跳過此記錄
        }

        int pos = record->pos + 1; // 將位置從 0-based 轉換為 1-based
        if (record->n_allele < 1) {
            std::cerr << "無效的等位基因數量\n";
            continue; // 如果等位基因數量無效，跳過此記錄
        }
        std::string ref_allele = record->d.allele[0]; // 獲取 REF 等位基因

        // 確保有 ALT 基因
        std::string alt_alleles;
        for (int i = 1; i < record->n_allele; i++) { // 可能有多個 ALT
            if (i > 1) alt_alleles += ","; // 如果有多個 ALT，使用逗號分隔
            alt_alleles += record->d.allele[i]; // 添加 ALT 等位基因
        }

        // 儲存到 `variants` 用於後續分析
        Variant var;
        var.chrom = chrom; // 設置染色體名稱
        var.pos = pos; // 設置位置
        var.ref = ref_allele; // 設置參考等位基因
        var.alt = alt_alleles; // 設置突變等位基因
        variants.push_back(var); // 將變異資訊添加到列表中
    }

    // 釋放資源
    bcf_destroy(record); // 釋放 VCF 讀取結構
    bcf_hdr_destroy(hdr); // 釋放標頭
    hts_close(vcfFile); // 關閉 VCF 文件
    return true; // 返回 true，表示成功載入變異
}

// 獲取變異資訊
const std::vector<Variant>& VCFParser::getVariants() const {
    return variants;
}