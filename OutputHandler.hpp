#ifndef OUTPUT_HANDLER_HPP
#define OUTPUT_HANDLER_HPP

#include <string>
#include <vector>
#include "Analysis.hpp"

class OutputHandler {
public:
    // 輸出 somatic 分析結果至 Somatic_analy.txt
    static bool writeSomaticAnaly(const std::vector<SomaticAnalyData> &somaticData,
                                  const std::string &outputFolder);
    // 輸出 CpG 甲基化分析結果至 methyl_analy.txt
    static bool writeMethylAnaly(const std::vector<MethylAnalyData> &methylData,
                                 const std::string &outputFolder);
};

#endif // OUTPUT_HANDLER_HPP
