#include <unordered_map>
#include <set>
#include <vector>
#include "BAMProcessor.hpp"

struct ReadInfoAggregated {
    std::string readName;
    std::set<int> mutatedPositions;      // 此 read 出現突變的所有絕對位置
    std::vector<int> allMethylationScores; // 聚合該 read 所有突變位點的甲基化分數
};

// 此函式接收所有變異位點分析得到的 read 資訊，並依 readName 整合
std::unordered_map<std::string, ReadInfoAggregated> aggregateReadInfo(const std::vector<std::vector<ReadInfo>> &allVariantReads) {
    std::unordered_map<std::string, ReadInfoAggregated> aggregated;
    for (const auto &variantReads : allVariantReads) {
        for (const auto &r : variantReads) {
            auto &agg = aggregated[r.readName];
            agg.readName = r.readName;
            // 若該 read 在此突變位點上帶突變，則記錄突變位置
            for (int pos : r.mutatedPositions) {
                agg.mutatedPositions.insert(pos);
            }
            // 累積甲基化分數
            agg.allMethylationScores.insert(agg.allMethylationScores.end(),
                                            r.methylationScores.begin(),
                                            r.methylationScores.end());
        }
    }
    return aggregated;
} 