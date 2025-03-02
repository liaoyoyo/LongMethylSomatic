#ifndef RESULTWRITER_HPP
#define RESULTWRITER_HPP

#include <string>
#include <vector>

class ResultWriter {
public:
    bool write(const std::string &outputFile, const std::vector<std::string>& results);
};

#endif 