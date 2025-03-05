#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include "CommonTypes.hpp"
#include <string>
#include <vector>

struct SomaticAnalyData {
    std::string chr;
    int pos;           // 1-based
    std::string ref;
    std::string alt;
    int ref_count;
    int alt_count;
    double ref_methyl;
    double alt_methyl;
};

struct MethylAnalyData {
    std::string chr;
    int pos;           // 1-based
    int somatic_pos;   // 1-based
    char somatic_base;
    double high_methyl;
};

struct AnalysisResult {
    std::vector<SomaticAnalyData> somaticData;
    std::vector<MethylAnalyData>  methylData;
};

class Analysis {
public:
    static AnalysisResult compute(const std::vector<SomaticSite>& somaticSites,
                                  const std::string &tumorBamFile,
                                  const std::string &normalBamFile,
                                  int window);
};

#endif
