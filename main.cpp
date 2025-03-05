#include "ArgParser.hpp"
#include "VCFHandler.hpp"
#include "Analysis.hpp"
#include "OutputHandler.hpp"
#include "Utility.hpp"
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[]) {
    // 解析命令列參數
    Args args = ArgParser::parse(argc, argv);
#ifdef _OPENMP
    omp_set_num_threads(args.maxThreads);
#endif

    // 總計時
    Timer totalTimer;

    // 解析 VCF 檔案，取得 somatic mutation 資訊
    std::cout << "開始解析 VCF 檔案..." << std::endl;
    Timer vcfTimer;
    auto somaticSites = VCFHandler::parseSomaticSites(args.vcfFile);
    std::cout << "VCF 解析耗時: " << vcfTimer.stop() << " 秒, 取得 " 
              << somaticSites.size() << " 筆 mutation" << std::endl;

    // 執行 BAM 分析：解析 CIGAR、MD 與 MM/ML 標籤，統計甲基化資訊
    std::cout << "開始執行分析..." << std::endl;
    Timer analysisTimer;
    AnalysisResult analysisResult = Analysis::compute(somaticSites, args.tumorBam, args.normalBam, args.window);
    std::cout << "分析耗時: " << analysisTimer.stop() << " 秒" << std::endl;

    // 輸出結果
    std::cout << "開始輸出結果..." << std::endl;
    Timer outputTimer;
    bool writeSomaticSuccess = OutputHandler::writeSomaticAnaly(analysisResult.somaticData, args.outputFolder);
    bool writeMethylSuccess = OutputHandler::writeMethylAnaly(analysisResult.methylData, args.outputFolder);
    if (!writeSomaticSuccess || !writeMethylSuccess) {
        std::cerr << "錯誤：結果輸出失敗" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "輸出耗時: " << outputTimer.stop() << " 秒" << std::endl;

    std::cout << "總執行時間: " << totalTimer.stop() << " 秒" << std::endl;
    return EXIT_SUCCESS;
}
