#ifndef BAMPROCESSOR_HPP
#define BAMPROCESSOR_HPP

#include <string>
#include <vector>
#include <htslib/sam.h>

// 用於存放 MM 標籤解析後的甲基化資訊
struct MethylationInfo {
    std::string mod_type;        // 例如 "C+m"、"C+h" 等
    std::vector<int> positions;  // 對應的 CpG 位點（以 read 中的相對位置或其他定義）
};

// 用於存放單條 read 的資訊
struct ReadInfo {
    std::string readName;                 // 讀取名稱
    bool isMutation;                      // 該 read 在當前突變位點上是否呈現突變（根據 ALT 比對）
    std::vector<MethylationInfo> methylationData; // 從 MM:Z: 標籤解析得到
    std::vector<int> methylationScores;   // 從 ML:B:C 標籤解析得到的分數
    std::vector<int> methylationPositions; // 解析出來的 CpG 位置（假設取自第一組 methylationData）
    std::vector<int> mutatedPositions;    // 紀錄該 read 上出現突變的所有突變位點（絕對位置）
};

class BAMProcessor {
public:
    BAMProcessor(const std::string &bamFilename);
    ~BAMProcessor();

    bool openBam();
    // 以指定染色體與突變位置、參考與突變碱基，並以 window 為範圍提取 read 資訊
    std::vector<ReadInfo> processVariant(const std::string &chr, int pos, char ref, char alt, int window);
    
    // 解析 MM:Z: 標籤，傳回各修飾類型與位置
    std::vector<MethylationInfo> parseMethylationTag(const std::string &mm_tag);
    // 解析 ML:B:C 標籤，傳回各分數
    std::vector<int> parseMLTag(uint8_t* mlTag);
    bam_hdr_t* getHeader() const;

private:
    std::string bamFilename;
    samFile* bamFile;
    bam_hdr_t* header;
    hts_idx_t* idx;
};

#endif // BAMPROCESSOR_HPP