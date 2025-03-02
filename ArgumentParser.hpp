#ifndef ARGUMENTPARSER_HPP
#define ARGUMENTPARSER_HPP

#include <string>
#include <iostream>

class ArgumentParser {
public:
    static ArgumentParser& getInstance();
    bool parse(int argc, char* argv[]);
    const std::string& getNormalBam() const;
    const std::string& getTumorBam() const;
    const std::string& getVcfFile() const;
    const std::string& getReference() const;
    const std::string& getOutputFile() const;
    int getThreadCount() const;
    void printUsage(const char* programName) const;
    std::string decompressGzFile(const std::string& gzFilePath);
private:
    ArgumentParser() = default;
    std::string normalBam;
    std::string tumorBam;
    std::string vcfFile;
    std::string reference;
    std::string outputFile;
    int threadCount = 4;
};

#endif 