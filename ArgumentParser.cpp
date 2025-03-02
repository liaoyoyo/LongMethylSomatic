#include "ArgumentParser.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <zlib.h>

// 獲取 ArgumentParser 的單例實例
ArgumentParser& ArgumentParser::getInstance() {
    static ArgumentParser instance; // 使用靜態變數確保單例
    return instance;
}

// 解析命令行參數
bool ArgumentParser::parse(int argc, char* argv[]) {
    if (argc < 11) { // 檢查參數數量
        printUsage(argv[0]);
        return false;
    }

    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-o" && i + 1 < argc) {
            outputFile = argv[++i]; // 獲取輸出文件名
        } else if (std::string(argv[i]) == "-n" && i + 1 < argc) {
            normalBam = argv[++i]; // 獲取正常 BAM 文件名
        } else if (std::string(argv[i]) == "-t" && i + 1 < argc) {
            tumorBam = argv[++i]; // 獲取腫瘤 BAM 文件名
        } else if (std::string(argv[i]) == "-v" && i + 1 < argc) {
            vcfFile = argv[++i]; // 獲取 VCF 文件名
        } else if (std::string(argv[i]) == "-r" && i + 1 < argc) {
            reference = argv[++i]; // 獲取參考基因組文件名
        } else if (std::string(argv[i]) == "-j" && i + 1 < argc) {
            threadCount = std::stoi(argv[++i]); // 獲取執行緒數量
        }
    }

    return true; // 參數解析成功
}

// 獲取正常 BAM 文件名
const std::string& ArgumentParser::getNormalBam() const {
    return normalBam;
}

// 獲取腫瘤 BAM 文件名
const std::string& ArgumentParser::getTumorBam() const {
    return tumorBam;
}

// 獲取 VCF 文件名
const std::string& ArgumentParser::getVcfFile() const {
    return vcfFile;
}

// 獲取參考基因組文件名
const std::string& ArgumentParser::getReference() const {
    return reference;
}

// 獲取輸出文件名
const std::string& ArgumentParser::getOutputFile() const {
    return outputFile;
}

// 獲取執行緒數量
int ArgumentParser::getThreadCount() const {
    return threadCount;
}

// 打印使用說明
void ArgumentParser::printUsage(const char* programName) const {
    std::cerr << "用法: " << programName 
              << " -o <output.txt> -n <normal.bam> -t <tumor.bam> -v <somatic.vcf> -r <reference.fa> -j <thread_count>\n";
}

// 解壓縮 .gz 文件
std::string ArgumentParser::decompressGzFile(const std::string& gzFilePath) {
    std::ifstream gzFile(gzFilePath, std::ios_base::in | std::ios_base::binary);
    if (!gzFile) {
        std::cerr << "無法打開文件: " << gzFilePath << "\n";
        return "";
    }

    std::ostringstream outStream;
    char buffer[128];
    gzFile.read(buffer, sizeof(buffer));
    while (gzFile) {
        outStream.write(buffer, gzFile.gcount());
        gzFile.read(buffer, sizeof(buffer));
    }

    return outStream.str(); // 返回解壓縮的內容
}