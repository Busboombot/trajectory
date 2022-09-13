
#include <chrono>
#include <thread>
#include <unistd.h>
#include "trj_util.h"
#include <iostream>
#include <string>
#include "col.h"

bool same_sign(float a, float b){
    return (a == 0) or (b == 0) or (sgn(a) == sgn(b));
}

int sign(int x) {
    if (x == 0) return 0;
    else if  (x > 0) return 1;
    else return -1;
}

int sign(trj_float_t x) {
    if (x == 0) return 0;
    else if  (x > 0) return 1;
    else return -1;
}


#ifdef TRJ_ENV_HOST
void delay(uint32_t ms){
    std::this_thread::sleep_for(std::chrono::milliseconds(ms));

}
void delayMicroseconds(uint32_t us){
    std::this_thread::sleep_for(std::chrono::microseconds(us));

}


steadyClock::time_point usince_start =  steadyClock::now();
void start_usince(){
    usince_start =  steadyClock::now();
}

uint32_t usince(){
    auto elapsed = steadyClock::now() - usince_start;
    return (uint32_t)(elapsed.count()/1000);
}
#else

uint32_t  usince_start =  micros();
void start_usince(){
     usince_start =  micros();
}

uint32_t usince(){
    return  micros() - usince_start;
}
#endif

std::vector<std::string> splitString(const std::string& str){
    std::vector<std::string> tokens;
 
    std::string::size_type pos = 0;
    std::string::size_type prev = 0;
    while ((pos = str.find('\n', prev)) != std::string::npos) {
        tokens.push_back(str.substr(prev, pos - prev));
        prev = pos + 1;
    }
    tokens.push_back(str.substr(prev));
 
    return tokens;
}


std::string yellow = col::make(col::yellow, col::def, false, false, false);
std::string green  = col::make(col::green,  col::def, false, false, false);
std::string blue  = col::make(col::blue,  col::def, false, false, false);
std::string blue_bg  = col::make(col::blue,  col::def, true, false, false);
std::string creset = "\x1b[0m";



