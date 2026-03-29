#include "TimeMeasure.hpp"

TimeMeasure::TimeMeasure(){
    _start = std::chrono::high_resolution_clock::now();
    _stop  = std::chrono::high_resolution_clock::now();
}

void TimeMeasure::start(){
    _start = std::chrono::high_resolution_clock::now();
    _stop  = std::chrono::high_resolution_clock::now();
}

void TimeMeasure::stop(){
    _stop = std::chrono::high_resolution_clock::now();
}

int TimeMeasure::sec(){
    auto e_delai = std::chrono::duration_cast<std::chrono::seconds>(_stop - _start);
    return  e_delai.count() == 0 ? 1 : e_delai.count();
}

int TimeMeasure::ms(){
    auto e_delai = std::chrono::duration_cast<std::chrono::milliseconds>(_stop - _start);
    return  e_delai.count() == 0 ? 1 : e_delai.count();
}

int TimeMeasure::us(){
    auto e_delai = std::chrono::duration_cast<std::chrono::microseconds>(_stop - _start);
    return  e_delai.count() == 0 ? 1 : e_delai.count();
}

int TimeMeasure::ns(){
    auto e_delai = std::chrono::duration_cast<std::chrono::nanoseconds>(_stop - _start);
    return  e_delai.count() == 0 ? 1 : e_delai.count();
}

int TimeMeasure::KBps(const int64_t bytes){
    const double t_ms   = ms();
    const double debit = bytes / t_ms;
    return (int)debit;
}

int TimeMeasure::MBps(const int64_t bytes){
    const double t_us   = us();
    const double debit = bytes / t_us;
    return (int)debit;
}

int TimeMeasure::GBps(const int64_t bytes){
    const double t_ns   = ns();
    const double debit = bytes / t_ns;
    return (int)debit;
}