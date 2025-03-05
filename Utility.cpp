#include "Utility.hpp"
#include <chrono>

Timer::Timer() {
    start();
}

void Timer::start() {
    start_time = std::chrono::high_resolution_clock::now();
}

double Timer::stop() {
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}
