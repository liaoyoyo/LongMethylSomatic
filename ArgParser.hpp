#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP

#include <string>

struct Args {
    std::string normalBam;    // -n 或 --normal (可選)
    std::string tumorBam;     // -t 或 --tumor (必填)
    std::string vcfFile;      // -v 或 --vcf (必填)
    std::string refFile;      // -r 或 --ref (可選)
    std::string outputFolder; // -o 或 --output (可選，預設為當前目錄)
    int window;               // -w 或 --window (可選，預設2000)
    int maxThreads;           // -j 或 --threads (可選，預設使用系統最大)
};

class ArgParser {
public:
    static Args parse(int argc, char* argv[]);
    static void printHelp(const char* progName);
    static void printArgs(const Args& args);
};

#endif // ARG_PARSER_HPP
