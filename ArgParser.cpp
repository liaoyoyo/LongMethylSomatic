#include "ArgParser.hpp"
#include <iostream>
#include <cstring>
#include <getopt.h>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif

// 顯示使用說明
void ArgParser::printHelp(const char* progName) {
    std::cout << "使用說明: " << progName << " [options]\n"
              << "選項:\n"
              << "  -n, --normal <file>    Normal BAM 檔案 (可選)\n"
              << "  -t, --tumor <file>     Tumor BAM 檔案 (必填)\n"
              << "  -v, --vcf <file>       Somatic VCF 檔案 (必填)\n"
              << "  -r, --ref <file>       參考基因組檔案 (可選)\n"
              << "  -o, --output <folder>  輸出資料夾 (預設 './')\n"
              << "  -w, --window <num>     Somatic 讀取範圍 (預設 2000)\n"
              << "  -j, --threads <num>    最大執行緒數 (預設使用最大值)\n"
              << "  -h, --help             顯示此訊息\n";
}

Args ArgParser::parse(int argc, char* argv[]) {
    Args args;
    // 設定預設值
    args.window = 2000;
    args.outputFolder = "./";
#ifdef _OPENMP
    args.maxThreads = omp_get_max_threads();
#else
    args.maxThreads = 1;
#endif

    static struct option longOptions[] = {
        {"normal", required_argument, 0, 'n'},
        {"tumor", required_argument, 0, 't'},
        {"vcf", required_argument, 0, 'v'},
        {"ref", required_argument, 0, 'r'},
        {"output", required_argument, 0, 'o'},
        {"window", required_argument, 0, 'w'},
        {"threads", required_argument, 0, 'j'},
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };

    int optionChar;
    int optionIndex = 0;
    while ((optionChar = getopt_long(argc, argv, "n:t:v:r:o:w:j:h", longOptions, &optionIndex)) != -1) {
        switch (optionChar) {
            case 'n':
                args.normalBam = optarg;
                break;
            case 't':
                args.tumorBam = optarg;
                break;
            case 'v':
                args.vcfFile = optarg;
                break;
            case 'r':
                args.refFile = optarg;
                break;
            case 'o':
                args.outputFolder = optarg;
                break;
            case 'w':
                args.window = std::stoi(optarg);
                break;
            case 'j':
                args.maxThreads = std::stoi(optarg);
                break;
            case 'h':
                printHelp(argv[0]);
                exit(EXIT_SUCCESS);
            case '?':
                printHelp(argv[0]);
                exit(EXIT_FAILURE);
            default:
                break;
        }
    }

    // 檢查必要參數
    if (args.tumorBam.empty()) {
        std::cerr << "錯誤：必須指定 Tumor BAM 檔案 (-t 或 --tumor)" << std::endl;
        printHelp(argv[0]);
        exit(EXIT_FAILURE);
    }
    if (args.vcfFile.empty()) {
        std::cerr << "錯誤：必須指定 Somatic VCF 檔案 (-v 或 --vcf)" << std::endl;
        printHelp(argv[0]);
        exit(EXIT_FAILURE);
    }
    return args;
}
