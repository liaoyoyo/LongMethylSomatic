#ifndef ANALYSISRESULT_HPP
#define ANALYSISRESULT_HPP

#include <string>

// 存儲突變位點與甲基化分析的結果
struct AnalysisResult {
    std::string somaticPosition;  // 例如 "chr1:1000 (A>G)"
    std::string cpgPosition;      // CpG 位置，例如 "chr1:950"
    int mutMethylated;            // 突變 reads 甲基化數
    int mutTotal;                 // 突變 reads 總數
    int nonMutMethylated;         // 非突變 reads 甲基化數
    int nonMutTotal;              // 非突變 reads 總數
    double score;                 // 加權分數
    std::string conclusion;        // 結論，例如 "存在一致性高甲基化差異"
};

#endif
