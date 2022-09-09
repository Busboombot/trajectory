
#include <chrono>
#include <thread>
#include <unistd.h>
#include "Arduino.h"
#include "trj_util.h"

bool same_sign(float a, float b){
    return (a == 0) or (b == 0) or (sgn(a) == sgn(b));
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