#include "ArgumentParser.hpp"
#include "VCFParser.hpp"
#include "BAMProcessor.hpp"
#include "MethylationAnalyzer.hpp"
#include "ResultWriter.hpp"
#include <iostream>
#include <vector>
#include <unordered_map>
#include <set>
#include <sstream>

// 前述 aggregateReadInfo 函式與相關結構
struct ReadInfoAggregated {
    std::string readName;
    std::set<int> mutatedPositions;
    std::vector<int> allMethylationScores;
};

std::unordered_map<std::string, ReadInfoAggregated> aggregateReadInfo(const std::vector<std::vector<ReadInfo>> &allVariantReads);

int main(int argc, char* argv[]) {
    // 解析命令列參數
    ArgumentParser &argParser = ArgumentParser::getInstance();
    if (!argParser.parse(argc, argv)) {
        return 1;
    }
    
    // 載入 VCF 檔案與變異資訊
    VCFParser vcfParser(argParser.getVcfFile());
    if (!vcfParser.loadVariants()) {
        std::cerr << "VCF 變異載入失敗，程式終止。\n";
        return 1;
    }
    const auto &variants = vcfParser.getVariants();
    if (variants.empty()) {
        std::cerr << "VCF 無突變資訊，程式終止。\n";
        return 1;
    }
    
    // 開啟腫瘤 BAM 檔案
    BAMProcessor tumorBamProcessor(argParser.getTumorBam());
    if (!tumorBamProcessor.openBam()) {
        return 1;
    }
    MethylationAnalyzer analyzer;
    
    // 儲存所有變異位點分析得到的 read 資訊
    std::vector<std::vector<ReadInfo>> allVariantReads;
    
    std::cout << "開始處理變異...\n";
    // 目標 1：針對每個突變位點，計算突變 reads 與非突變 reads 的甲基化比例
    for (const auto &var : variants) {
        auto reads = tumorBamProcessor.processVariant(var.chrom, var.pos, var.ref[0], var.alt[0], 50);
        if (reads.empty()) {
            std::cerr << "未找到任何 reads 於變異: " << var.chrom << ":" << var.pos << "\n";
            continue; // 跳過此變異
        }
        allVariantReads.push_back(reads);
        
        double mutMethRate = analyzer.analyze(reads, true);
        double nonMutMethRate = analyzer.analyze(reads, false);
        
        std::cout << "處理變異: " << var.chrom << ":" << var.pos << "\n";
        std::cout << var.chrom << ":" << var.pos << " (" << var.ref << ">" << var.alt << ")\n";
        std::cout << "突變 reads 甲基化比例: " << mutMethRate << "%\n";
        std::cout << "非突變 reads 甲基化比例: " << nonMutMethRate << "%\n";
    }
    std::cout << "變異處理完成，開始聚合 reads...\n";
    
    // 目標 2：彙整同一 read 在不同突變位點的資訊
    auto aggregatedReads = aggregateReadInfo(allVariantReads);
    std::cout << "\nRead-level 統計資訊：\n";
    std::vector<std::string> resultStrings;
    for (const auto &entry : aggregatedReads) {
        const auto &agg = entry.second;
        double avgMeth = 0.0;
        if (!agg.allMethylationScores.empty()) {
            double sum = 0;
            for (auto score : agg.allMethylationScores) sum += score;
            avgMeth = sum / agg.allMethylationScores.size();
        }
        std::ostringstream oss;
        oss << "Read: " << agg.readName
            << ", 突變位點數: " << agg.mutatedPositions.size()
            << ", 平均甲基化分數: " << avgMeth;
        resultStrings.push_back(oss.str());
    }
    
    // 在計算完所有變異後，將結果寫入文件
    ResultWriter writer;
    if (!writer.write(argParser.getOutputFile(), resultStrings)) {
        std::cerr << "結果輸出失敗\n";
        return 1;
    }
    
    return 0;
} 