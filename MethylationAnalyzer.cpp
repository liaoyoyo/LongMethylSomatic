#include "MethylationAnalyzer.hpp"
#include <iostream>

double MethylationAnalyzer::analyze(const std::vector<ReadInfo>& reads, bool isMutation) {
    int totalReads = 0;
    double totalScore = 0.0;

    // 計算突變reads數量
    for (const auto &r : reads) {
        if (r.isMutation == isMutation) {
            ++totalReads;
            // 假設對整個 read 的甲基化分數做平均
            double sum = 0.0;
            for (int val : r.methylationScores) {
                sum += val;
            }
            if (!r.methylationScores.empty()) {
                double avg = sum / r.methylationScores.size();
                totalScore += avg; 
            }
        }
    }

    if (totalReads > 0) {
        return totalScore / totalReads; // 返回平均分數
    } else {
        std::cerr << "無讀取數據，返回 0\n"; // 添加調試信息
        return 0.0; // 如果沒有讀取，返回 0
    }
} 