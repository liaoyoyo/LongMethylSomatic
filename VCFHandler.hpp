#ifndef VCF_HANDLER_HPP
#define VCF_HANDLER_HPP

#include "CommonTypes.hpp"
#include <string>
#include <vector>

class VCFHandler {
public:
    // 解析 VCF 檔案，回傳所有 somatic mutation 位點資訊
    static std::vector<SomaticSite> parseSomaticSites(const std::string &vcfFile);
};

#endif // VCF_HANDLER_HPP
