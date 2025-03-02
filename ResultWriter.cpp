#include "ResultWriter.hpp"
#include <fstream>
#include <iostream>

// 寫入結果到指定的輸出文件
bool ResultWriter::write(const std::string &outputFile, const std::vector<std::string>& results) {
    std::ofstream ofs(outputFile);
    if (!ofs) {
        std::cerr << "錯誤: 無法開啟輸出檔案 " << outputFile << ".\n";
        return false; // 返回失敗
    }
    for (const auto &line : results) {
        ofs << line << "\n"; // 寫入每一行結果
    }
    ofs.close(); // 關閉文件
    return true; // 返回成功
}