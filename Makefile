CXX = g++
CXXFLAGS = -std=c++11 -g -O2 -Wall -fopenmp -I/big8_disk/liaoyoyo2001/test1/htslib/include
LDFLAGS = -L/big8_disk/liaoyoyo2001/test1/htslib/lib
LIBS = -lhts -llzma -lz -lbz2 -ldeflate -lcurl -lssl -lcrypto -lpthread -lm 

SRCS = main.cpp ArgParser.cpp VCFHandler.cpp Analysis.cpp OutputHandler.cpp Utility.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = LongMethylSomatic

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
