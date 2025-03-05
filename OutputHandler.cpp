#include "OutputHandler.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>

// 建立輸出資料夾（跨平台）
static void createDirectory(const std::string &folder) {
#ifdef _WIN32
    _mkdir(folder.c_str());
#else
    mkdir(folder.c_str(), 0755);
#endif
}

// 輸出 somatic mutation 分析結果
bool OutputHandler::writeSomaticAnaly(const std::vector<SomaticAnalyData> &somaticData,
                                      const std::string &outputFolder) {
    createDirectory(outputFolder);
    std::string filename = outputFolder + "/Somatic_analy.txt";
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "錯誤：無法開啟檔案 " << filename << " 進行輸出" << std::endl;
        return false;
    }
    ofs << "chr\tPOS\tref\talt\tref_count\talt_count\tref_methyl\talt_methyl\n";
    for (const auto &data : somaticData) {
        ofs << data.chr << "\t" << data.pos << "\t" << data.ref << "\t"
            << data.alt << "\t" << data.ref_count << "\t"
            << data.alt_count << "\t" << data.ref_methyl << "\t"
            << data.alt_methyl << "\n";
    }
    ofs.close();
    std::cout << "Somatic_analy.txt 輸出完成" << std::endl;
    return true;
}

// 輸出 methylation 分析結果
bool OutputHandler::writeMethylAnaly(const std::vector<MethylAnalyData> &methylData,
                                     const std::string &outputFolder) {
    // 依染色體、位點、somatic 位點與 observed allele排序；methylation score 由高到低排列
    std::vector<MethylAnalyData> sortedData = methylData;
    std::sort(sortedData.begin(), sortedData.end(), [](const MethylAnalyData &a, const MethylAnalyData &b) {
        if (a.chr != b.chr)
            return a.chr < b.chr;
        if (a.pos != b.pos)
            return a.pos < b.pos;
        if (a.somatic_pos != b.somatic_pos)
            return a.somatic_pos < b.somatic_pos;
        if (a.somatic_base != b.somatic_base)
            return a.somatic_base < b.somatic_base;
        return a.high_methyl > b.high_methyl;
    });

    createDirectory(outputFolder);
    std::string filename = outputFolder + "/methyl_analy.txt";
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "錯誤：無法開啟檔案 " << filename << " 進行輸出" << std::endl;
        return false;
    }
    ofs << "Methyl_Chr\tMethyl_POS\tSomatic_POS\tSomatic_Allele\tMethylation_Score\n";
    for (const auto &data : sortedData) {
        ofs << data.chr << "\t" << data.pos << "\t"
            << data.somatic_pos << "\t" << data.somatic_base << "\t"
            << data.high_methyl << "\n";
    }
    ofs.close();
    std::cout << "methyl_analy.txt 輸出完成" << std::endl;
    return true;
}
