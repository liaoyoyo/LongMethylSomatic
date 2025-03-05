#include "VCFHandler.hpp"
#include <iostream>
#include "htslib/vcf.h"
#include "htslib/hts.h"

std::vector<SomaticSite> VCFHandler::parseSomaticSites(const std::string &vcfFile) {
    std::vector<SomaticSite> sites;
    
    // 開啟 VCF 檔案
    htsFile *vcfFp = bcf_open(vcfFile.c_str(), "r");
    if (!vcfFp) {
        std::cerr << "錯誤：無法開啟 VCF 檔案 " << vcfFile << std::endl;
        exit(EXIT_FAILURE);
    }
    bcf_hdr_t *header = bcf_hdr_read(vcfFp);
    bcf1_t *record = bcf_init();
    
    // 逐筆讀取 VCF 記錄
    while (bcf_read(vcfFp, header, record) == 0) {
        bcf_unpack(record, BCF_UN_ALL);
        SomaticSite site;
        site.chr = bcf_hdr_id2name(header, record->rid);
        site.pos = record->pos + 1; // VCF pos 為 0-based，轉為 1-based
        site.ref = record->d.allele[0];
        site.alt = (record->d.m_allele > 1) ? record->d.allele[1] : "";
        sites.push_back(site);
    }
    
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    bcf_close(vcfFp);
    return sites;
}
