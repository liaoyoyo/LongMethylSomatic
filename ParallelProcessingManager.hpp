#ifndef PARALLELPROCESSINGMANAGER_HPP
#define PARALLELPROCESSINGMANAGER_HPP

#include "VCFParser.hpp"
#include "BAMProcessor.hpp"
#include "MethylationAnalyzer.hpp"
#include "ResultWriter.hpp"
#include "AnalysisResult.hpp"
#include <vector>
#include <string>

class ParallelProcessingManager {
public:
    // 構造函數，初始化變異列表和處理器
    ParallelProcessingManager(const std::vector<Variant>& vars,
                              BAMProcessor* normalBamProc,
                              BAMProcessor* tumorBamProc,
                              MethylationAnalyzer* methAnalyzer);
    void process(); // 處理變異
    const std::vector<AnalysisResult>& getResults() const; // 獲取分析結果
    std::vector<std::string> getResultsAsStrings() const; // 獲取格式化的結果字符串

    // 全域變數
    int windowSize; // 窗口大小
    int cpgOffset;  // CpG 偏移量

private:
    std::vector<Variant> variants; // 存儲變異列表
    BAMProcessor* normalBamProcessor; // 正常 BAM 處理器
    BAMProcessor* tumorBamProcessor; // 腫瘤 BAM 處理器
    MethylationAnalyzer* methylationAnalyzer; // 甲基化分析器
    std::vector<AnalysisResult> results; // 存儲分析結果
};

#endif 