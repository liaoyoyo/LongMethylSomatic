#include "Analysis.hpp"
#include "htslib/sam.h"  // 用於解析 BAM 標籤
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <vector>
#include <sstream>
#include <string>
#include <tuple>
#include <map>

//--------------------------------------------------
// 將 read 序列對應到參考座標 (1-based)
//--------------------------------------------------
static std::vector<int> buildReadToRefMap(const bam1_t* aln) {
    int readLength = aln->core.l_qseq;
    std::vector<int> readToRef(readLength + 1, -1);
    int refPos = aln->core.pos + 1; // 1-based
    uint32_t *cigar = bam_get_cigar(aln);
    int nCigar = aln->core.n_cigar;
    int readPos = 1; // 1-based read 座標

    for (int i = 0; i < nCigar; i++) {
        int op = cigar[i] & BAM_CIGAR_MASK;
        int len = cigar[i] >> BAM_CIGAR_SHIFT;
        switch (op) {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                for (int j = 0; j < len; j++) {
                    if (readPos <= readLength)
                        readToRef[readPos] = refPos;
                    readPos++;
                    refPos++;
                }
                break;
            case BAM_CINS:
            case BAM_CSOFT_CLIP:
                readPos += len;
                break;
            case BAM_CDEL:
            case BAM_CREF_SKIP:
                refPos += len;
                break;
            default:
                break;
        }
    }
    return readToRef;
}

//--------------------------------------------------
// 解析 BAM 中的 MM/ML 標籤，取得甲基化記錄
//--------------------------------------------------
struct MethylationRecord {
    int refPos;    // 1-based 參考座標
    double prob;   // 修飾機率 (0 ~ 1)
};

static std::vector<MethylationRecord> parseMethylation(const bam1_t* aln) {
    std::vector<MethylationRecord> records;
    hts_base_mod_state *modState = hts_base_mod_state_alloc();
    int ret = bam_parse_basemod(aln, modState);
    if (ret < 0) {
        hts_base_mod_state_free(modState);
        return records;
    }
    auto readToRef = buildReadToRefMap(aln);
    int readLength = aln->core.l_qseq;
    hts_base_mod mods[5];
    int readPos0;
    int n = bam_next_basemod(aln, modState, mods, 5, &readPos0);
    while (n > 0) {
        for (int i = 0; i < n; i++) {
            int readPos1 = readPos0 + 1; // 轉為 1-based
            if (readPos1 >= 1 && readPos1 <= readLength) {
                int refPos = readToRef[readPos1];
                if (refPos > 0) {
                    double prob = mods[i].qual / 255.0;
                    records.push_back({refPos, prob});
                }
            }
        }
        n = bam_next_basemod(aln, modState, mods, 5, &readPos0);
    }
    hts_base_mod_state_free(modState);
    return records;
}

//--------------------------------------------------
// 分析每個 somatic site，統計 allele 與甲基化資訊
//--------------------------------------------------
AnalysisResult Analysis::compute(const std::vector<SomaticSite>& somaticSites,
                                   const std::string &tumorBamFile,
                                   const std::string &normalBamFile,
                                   int window) {
    AnalysisResult result;
    result.somaticData.resize(somaticSites.size());
    
    // 使用 OpenMP 平行處理各個 somatic site
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(somaticSites.size()); i++) {
        const auto& site = somaticSites[i];
        // 每個 thread 分別開啟 BAM 檔案
        samFile *bamFile = sam_open(tumorBamFile.c_str(), "r");
        if (!bamFile) continue;
        bam_hdr_t *header = sam_hdr_read(bamFile);
        hts_idx_t *index = sam_index_load(bamFile, tumorBamFile.c_str());
        if (!index) {
            bam_hdr_destroy(header);
           	sam_close(bamFile);
            continue;
        }
        // 定義查詢區間，確保範圍不小於 1
        int regionStart = (site.pos - window > 0) ? site.pos - window : 1;
        int regionEnd = site.pos + window;
        std::string regionStr = site.chr + ":" + std::to_string(regionStart) + "-" + std::to_string(regionEnd);
        hts_itr_t *iter = sam_itr_querys(index, header, regionStr.c_str());
        if (!iter) {
            hts_idx_destroy(index);
            bam_hdr_destroy(header);
            sam_close(bamFile);
            continue;
        }
        
        int refCount = 0, altCount = 0;
        double refMethSum = 0.0, altMethSum = 0.0;
        std::vector<MethylAnalyData> localMethyl; // 儲存此區域的 methylation 資料
        
        bam1_t *aln = bam_init1();
        while (sam_itr_next(bamFile, iter, aln) >= 0) {
            // 僅處理主對齊
            if ((aln->core.flag & BAM_FSECONDARY) || (aln->core.flag & BAM_FSUPPLEMENTARY))
                continue;
                
            // 取得此 read 的甲基化記錄
            auto methRecords = parseMethylation(aln);
            std::vector<MethylationRecord> filteredRecords;
            // 過濾僅保留在 somatic site ± window 範圍內的記錄
            for (const auto& rec : methRecords) {
                if (std::abs(rec.refPos - site.pos) <= window)
                    filteredRecords.push_back(rec);
            }
            double avgMethyl = 0.0;
            if (!filteredRecords.empty()) {
                double sumProb = 0.0;
                for (const auto& rec : filteredRecords)
                    sumProb += rec.prob;
                avgMethyl = sumProb / filteredRecords.size();
            }
            
            // 取得 read 在 somatic site 上的對應位置
            int readLength = aln->core.l_qseq;
            auto readToRef = buildReadToRefMap(aln);
            int readPosAtSite = -1;
            for (int pos = 1; pos <= readLength; pos++) {
                if (readToRef[pos] == site.pos) {
                    readPosAtSite = pos;
                    break;
                }
            }
            char observedAllele = 'N';
            if (readPosAtSite > 0) {
                uint8_t *seq = bam_get_seq(aln);
                int baseVal = bam_seqi(seq, readPosAtSite - 1);
                observedAllele = "=ACMGRSVTWYHKDBN"[baseVal];
            }
            
            // 統計 ref 與 alt 的讀數及甲基化累計
            if (toupper(observedAllele) != toupper(site.ref[0])) {
                altCount++;
                altMethSum += avgMethyl;
            } else {
                refCount++;
                refMethSum += avgMethyl;
            }
            
            // 記錄所有符合條件的 methylation 資料
            for (const auto& rec : filteredRecords) {
                MethylAnalyData md;
                md.chr = site.chr;
                md.pos = rec.refPos;
                md.somatic_pos = site.pos;
                md.somatic_base = observedAllele;
                md.high_methyl = rec.prob;
                localMethyl.push_back(md);
            }
        }
        bam_destroy1(aln);
        hts_itr_destroy(iter);
        hts_idx_destroy(index);
        bam_hdr_destroy(header);
        sam_close(bamFile);
        
        double refMethyl = (refCount > 0) ? refMethSum / refCount : 0.0;
        double altMethyl = (altCount > 0) ? altMethSum / altCount : 0.0;
        
        // 將結果寫入全域結果，使用 critical 區段避免多緒競爭
        #pragma omp critical
        {
            result.somaticData[i] = {site.chr, site.pos, site.ref, site.alt, refCount, altCount, refMethyl, altMethyl};
            result.methylData.insert(result.methylData.end(), localMethyl.begin(), localMethyl.end());
        }
    } // end parallel for

    // 聚合相同位點的 methylation 資料
    {
        using MethylKey = std::tuple<std::string, int, int, char>;
        std::map<MethylKey, std::pair<double, int>> aggMap;
        for (const auto &m : result.methylData) {
            MethylKey key = std::make_tuple(m.chr, m.pos, m.somatic_pos, m.somatic_base);
            aggMap[key].first += m.high_methyl;
            aggMap[key].second += 1;
        }
        std::vector<MethylAnalyData> aggregated;
        for (const auto &entry : aggMap) {
            MethylAnalyData md;
            md.chr = std::get<0>(entry.first);
            md.pos = std::get<1>(entry.first);
            md.somatic_pos = std::get<2>(entry.first);
            md.somatic_base = std::get<3>(entry.first);
            md.high_methyl = entry.second.first / entry.second.second;
            aggregated.push_back(md);
        }
        result.methylData = std::move(aggregated);
    }
    
    if (!normalBamFile.empty()) {
        #pragma omp critical
        {
            std::cerr << "注意：Normal BAM 分析尚未實作。\n";
        }
    }
    return result;
}
