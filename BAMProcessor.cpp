#include "BAMProcessor.hpp"
#include <iostream>
#include <sstream>
#include <cstdlib>

// 建構子與解構子
BAMProcessor::BAMProcessor(const std::string &bamFilename)
    : bamFilename(bamFilename), bamFile(nullptr), header(nullptr), idx(nullptr) {}

BAMProcessor::~BAMProcessor() {
    if (idx) hts_idx_destroy(idx);
    if (header) bam_hdr_destroy(header);
    if (bamFile) sam_close(bamFile);
}

bool BAMProcessor::openBam() {
    bamFile = sam_open(bamFilename.c_str(), "r");
    if (!bamFile) {
        std::cerr << "無法開啟 BAM 檔案: " << bamFilename << "\n";
        return false;
    }
    header = sam_hdr_read(bamFile);
    if (!header) {
        std::cerr << "無法讀取 BAM header: " << bamFilename << "\n";
        sam_close(bamFile);
        return false;
    }
    idx = sam_index_load(bamFile, bamFilename.c_str());
    if (!idx) {
        std::cerr << "無法載入 BAM 索引: " << bamFilename << "\n";
        sam_close(bamFile);
        return false;
    }
    return true;
}

std::vector<MethylationInfo> BAMProcessor::parseMethylationTag(const std::string &mm_tag) {
    std::vector<MethylationInfo> methylation_data;
    std::stringstream ss(mm_tag);
    std::string token;
    while (std::getline(ss, token, ';')) {
        if (token.empty()) continue;
        std::stringstream sub_ss(token);
        std::string mod_type;
        std::getline(sub_ss, mod_type, ',');
        MethylationInfo meth_info;
        meth_info.mod_type = mod_type;
        std::string pos_str;
        while (std::getline(sub_ss, pos_str, ',')) {
            try {
                meth_info.positions.push_back(std::stoi(pos_str));
            } catch (...) {
                std::cerr << "解析 MM 標籤錯誤: " << pos_str << "\n";
            }
        }
        methylation_data.push_back(meth_info);
    }
    return methylation_data;
}

std::vector<int> BAMProcessor::parseMLTag(uint8_t* mlTag) {
    std::vector<int> likelihoods;
    if (!mlTag) return likelihoods;
    uint8_t* ptr = mlTag;

    // 確認格式：應以 'B' 與 'C' 開頭
    if (*ptr != 'B' || *(ptr + 1) != 'C') {
        std::cerr << "ML 標籤格式錯誤\n";
        return likelihoods;
    }

    // 根據實際格式跳過前置資訊，此處假設跳過 5 個位元組（包含 "B", "C" 與長度資訊）
    ptr += 5;

    // 取得數據長度（假設存為 32 位元 little-endian，從 mlTag+2 取得）
    int length = *((int32_t*)(mlTag + 2));

    for (int i = 0; i < length; ++i) {
        likelihoods.push_back(ptr[i]);
    }
    return likelihoods;
}

std::vector<ReadInfo> BAMProcessor::processVariant(const std::string &chr, int pos, char ref, char alt, int window) {
    std::vector<ReadInfo> readInfos;
    int tid = bam_name2id(header, chr.c_str());
    if (tid < 0) {
        std::cerr << "找不到染色體: " << chr << "\n";
        return readInfos;
    }

    // 確保變異位置在 BAM 的長度範圍內
    if (pos < 0 || pos >= header->target_len[tid]) {
        std::cerr << "變異位置超出 BAM 長度: " << chr << ":" << pos << "\n";
        return readInfos;
    }

    int beg = (pos - window - 1) < 0 ? 0 : pos - window - 1;
    int end = (pos + window) > header->target_len[tid] ? header->target_len[tid] : pos + window;

    hts_itr_t* iter = sam_itr_queryi(idx, tid, beg, end);
    bam1_t* aln = bam_init1();
    while (sam_itr_next(bamFile, iter, aln) >= 0) {
        ReadInfo info;
        // 取得 read 名稱
        info.readName = bam_get_qname(aln);
        
        // 解析 MM 與 ML 標籤取得甲基化資訊
        uint8_t* mmTag = bam_aux_get(aln, "MM");
        if (mmTag) {
            std::string mm_str = bam_aux2Z(mmTag);
            info.methylationData = parseMethylationTag(mm_str);
            // 假設取第一組 methylationData 作為 CpG 位點
            if (!info.methylationData.empty())
                info.methylationPositions = info.methylationData[0].positions;
        }
        uint8_t* mlTag = bam_aux_get(aln, "ML");
        if (mlTag) {
            info.methylationScores = parseMLTag(mlTag);
        }
        
        // 判斷該 read 在突變位點的狀態
        int read_offset = pos - aln->core.pos;  // 注意：aln->core.pos 為 0-based，pos 為 1-based
        if (read_offset >= 0 && read_offset < aln->core.l_qseq) {
            uint8_t* seq = bam_get_seq(aln);
            int base = bam_seqi(seq, read_offset);
            char base_char = seq_nt16_str[base];
            info.isMutation = (base_char == alt);
            if (info.isMutation) {
                // 記錄該突變的絕對位置
                info.mutatedPositions.push_back(pos);
            }
        } else {
            info.isMutation = false;
        }
        readInfos.push_back(info);
    }
    bam_destroy1(aln);
    hts_itr_destroy(iter);
    return readInfos;
}

bam_hdr_t* BAMProcessor::getHeader() const {
    return header; // 返回 header
}