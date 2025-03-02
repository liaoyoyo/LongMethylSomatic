# Makefile for VCFParser and somatic_methylation_analyzer

# 編譯器
CXX = g++
# 編譯選項
CXXFLAGS = -std=c++11 -fopenmp -g
# 連結庫
LIBS = -lhts
# 目標可執行檔
TARGET = somatic_methylation_analyzer
# 所有的源碼檔案
SOURCES = $(wildcard *.cpp)
# 生成的物件檔案
OBJECTS = $(SOURCES:.cpp=.o)

# 預設目標
all: $(TARGET)

# 編譯可執行檔
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# 編譯源碼檔案
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 清理生成的檔案
clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean 