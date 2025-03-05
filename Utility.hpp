#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <chrono>

class Timer {
public:
    Timer();
    void start();
    double stop();
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
};

#endif // UTILITY_HPP
