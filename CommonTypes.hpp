#ifndef COMMON_TYPES_HPP
#define COMMON_TYPES_HPP

#include <string>
#include <vector>

// 用來儲存 somatic mutation 資訊
struct SomaticSite {
    std::string chr;
    int pos;         // 1-based
    std::string ref;
    std::string alt;
};

#endif // COMMON_TYPES_HPP 