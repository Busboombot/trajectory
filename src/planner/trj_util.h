#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdint> 
#include <math.h> 
#include <chrono>
#include <stdexcept>

typedef std::chrono::milliseconds ms;
typedef std::chrono::microseconds us;
typedef std::chrono::steady_clock steadyClock;
typedef std::chrono::duration<uint64_t, std::micro> duration;

long micros();
void start_usince();
uint32_t usince();

#ifdef TRJ_ENV_HOST
void delay(uint32_t ms);
void delayMicroseconds(uint32_t us);
#endif

extern int here_count;
#define HERE(x) std::cout << "!!!HERE!!! " << x << " " << here_count++ << std::endl;

// Convert seconds to ticks, usually microseconds
#define SEC_TO_TICKS(v) (static_cast<int>(rint(v*1e6)))

using IntVec = std::vector<int32_t> ;

// Return -1 , 0  or 1 to indicate the sign of the number. 
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

bool same_sign(float a, float b);
int sign(int a);

// To convert class enums to ints. 
//template <typename E>
//constexpr typename std::underlying_type<E>::type to_underlying(E e) noexcept {
//    return static_cast<typename std::underlying_type<E>::type>(e);
//};

std::vector<std::string> splitString(const std::string& str);

void hello();
