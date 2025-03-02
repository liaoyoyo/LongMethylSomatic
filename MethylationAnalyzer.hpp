#ifndef METHYLATIONANALYZER_HPP
#define METHYLATIONANALYZER_HPP

#include "BAMProcessor.hpp"
#include <vector>

class MethylationAnalyzer {
public:
    // 分析提供的 reads 並回傳甲基化比例 (此處僅為模擬邏輯)
    double analyze(const std::vector<ReadInfo>& reads, bool isMutation);
};

#endif 