#pragma once
//
#include <cstdint>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <string>
#include <chrono>

class TimeMeasure{
private:
    std::chrono::time_point<std::chrono::steady_clock> _start;
    std::chrono::time_point<std::chrono::steady_clock> _stop;

public:
    TimeMeasure();

    void start();
    void stop();

    int sec();
    int ms();
    int us();
    int ns();

    int KBps(const int64_t bytes);
    int MBps(const int64_t bytes);
    int GBps(const int64_t bytes);
};
