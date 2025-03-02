#include "ParallelProcessingManager.hpp"
#include "VCFParser.hpp"
#include "BAMProcessor.hpp"
#include "MethylationAnalyzer.hpp"
#include "ResultWriter.hpp"
#include "AnalysisResult.hpp"
#include "ArgumentParser.hpp"
#include <vector>
#include <string>
#include <omp.h>
#include <sstream>
#include <iostream>
#include <algorithm>

ParallelProcessingManager::ParallelProcessingManager(const std::vector<Variant>& vars,
                                                     BAMProcessor* normalBamProc,
                                                     BAMProcessor* tumorBamProc,
                                                     MethylationAnalyzer* methAnalyzer)
    : variants(vars), normalBamProcessor(normalBamProc), tumorBamProcessor(tumorBamProc), methylationAnalyzer(methAnalyzer) {
    windowSize = 100000; // 設定窗口大小
    cpgOffset = 100000;  // 設定 CpG 偏移量
}

void ParallelProcessingManager::process() {
    results.resize(variants.size());
    int threadCount = ArgumentParser::getInstance().getThreadCount(); // 獲取執行緒數量
    #pragma omp parallel for schedule(dynamic) num_threads(threadCount) // 使用指定的執行緒數量
    for (size_t i = 0; i < variants.size(); i++) {
        const Variant &var = variants[i];

        int adjustedPos = var.pos; // 預設為原始位置

        // 確保變異位置不小於 0
        if (adjustedPos < 0) {
            std::cerr << "變異位置無效: " << var.chrom << ":" << adjustedPos << "\n";
            continue; // 跳過此變異
        }

        // 確保變異位置不大於 BAM 的長度
        if (adjustedPos > tumorBamProcessor->getHeader()->target_len[bam_name2id(tumorBamProcessor->getHeader(), var.chrom.c_str())]) {
            std::cerr << "變異位置超出 BAM 長度: " << var.chrom << ":" << adjustedPos << "\n";
            continue; // 跳過此變異
        }

        // 處理突變 reads
        std::vector<ReadInfo> tumorReads = tumorBamProcessor->processVariant(var.chrom, var.pos, var.ref[0], var.alt[0], windowSize);

        // 計算甲基化比例
        double mutMeth = methylationAnalyzer->analyze(tumorReads, true);
        // 計算非突變甲基化比例
        double nonMutMeth = methylationAnalyzer->analyze(tumorReads, false);
        // 計算突變 reads 數量
        int mutTotal = std::count_if(tumorReads.begin(), tumorReads.end(), [](const ReadInfo& r) { return r.isMutation; });
        // 計算非突變 reads 數量
        int nonMutTotal = tumorReads.size() - mutTotal;

        AnalysisResult res;
        res.somaticPosition = var.chrom + ":" + std::to_string(adjustedPos) + " (" + var.ref + ">" + var.alt + ")";
        res.cpgPosition = var.chrom + ":" + std::to_string(adjustedPos - cpgOffset); // 使用全域變數 cpgOffset
        res.mutMethylated = static_cast<int>(mutMeth * mutTotal / 100);
        res.mutTotal = mutTotal;
        res.nonMutMethylated = static_cast<int>(nonMutMeth * nonMutTotal / 100);
        res.nonMutTotal = nonMutTotal;
        res.score = std::abs(mutMeth - nonMutMeth);
        res.conclusion = (res.score > 50) ? "存在一致性高甲基化差異" : "無顯著一致性差異";
        results[i] = res;
    }
}

std::vector<std::string> ParallelProcessingManager::getResultsAsStrings() const {
    std::vector<std::string> output;
    for (const auto &res : results) {
        std::ostringstream oss;
        oss << res.somaticPosition << "\t" 
            << res.cpgPosition << "\t"
            << "突變 reads 甲基化比例: " << (res.mutTotal > 0 ? (100.0 * res.mutMethylated / res.mutTotal) : 0) << "% ("
            << res.mutMethylated << "/" << res.mutTotal << ")\t"
            << "非突變 reads 甲基化比例: " << (res.nonMutTotal > 0 ? (100.0 * res.nonMutMethylated / res.nonMutTotal) : 0) << "% ("
            << res.nonMutMethylated << "/" << res.nonMutTotal << ")\t"
            << "加權分數: " << res.score << "\t"
            << "結論: " << res.conclusion;
        output.push_back(oss.str());
    }
    return output;
}

const std::vector<AnalysisResult>& ParallelProcessingManager::getResults() const {
    return results;
} 