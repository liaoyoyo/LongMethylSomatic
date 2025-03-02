#ifndef VCFPARSER_HPP
#define VCFPARSER_HPP

#include <vector>
#include <string>
#include <htslib/vcf.h> // htslib VCF 解析
#include <htslib/hts.h> // htslib 文件處理

// 儲存突變資訊結構
struct Variant {
    std::string chrom;
    int pos; // 1-based 位置
    std::string ref;
    std::string alt;
};

class VCFParser {
public:
    explicit VCFParser(const std::string &filename);
    bool loadVariants();
    const std::vector<Variant>& getVariants() const;
private:
    std::string filename;
    std::vector<Variant> variants;
};

#endif
