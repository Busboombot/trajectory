#pragma once
#include <iostream>



//#define DEBUG_PIN_1 5
//#define DEBUG_PIN_2 6
//#define DEBUG_PIN_3 7
//#define DEBUG_PIN_4 8

#ifdef TRJ_ENV_HOST
#define DEBUG_SET_1
#define DEBUG_CLEAR_1
#define DEBUG_TOG_1 
#define DEBUG_SET_2
#define DEBUG_CLEAR_2
#define DEBUG_TOG_2
#define DEBUG_SET_3
#define DEBUG_CLEAR_3
#define DEBUG_TOG_3
#define DEBUG_SET_4
#define DEBUG_CLEAR_4
#define DEBUG_TOG_4

#else //TRJ_ENV_HOST

#ifdef DEBUG_PIN_1
extern int debug_state_1;
#define DEBUG_SET_1 digitalWriteFast(DEBUG_PIN_1, HIGH);
#define DEBUG_CLEAR_1 digitalWriteFast(DEBUG_PIN_1, LOW);
#define DEBUG_TOG_1 digitalWriteFast(DEBUG_PIN_1, debug_state_1=!debug_state_1);
#else
#define DEBUG_SET_1
#define DEBUG_CLEAR_1
#define DEBUG_TOG_1 
#endif

#ifdef DEBUG_PIN_2
extern int debug_state_2;
#define DEBUG_SET_2 digitalWriteFast(DEBUG_PIN_2, HIGH);
#define DEBUG_CLEAR_2 digitalWriteFast(DEBUG_PIN_2, LOW);
#define DEBUG_TOG_2 digitalWriteFast(DEBUG_PIN_2, debug_state_2=!debug_state_2);
#else
#define DEBUG_SET_2
#define DEBUG_CLEAR_2
#define DEBUG_TOG_2
#endif

#ifdef DEBUG_PIN_3
extern int debug_state_3;
#define DEBUG_SET_3 digitalWriteFast(DEBUG_PIN_3, HIGH);
#define DEBUG_CLEAR_3 digitalWriteFast(DEBUG_PIN_3, LOW);
#define DEBUG_TOG_3 digitalWriteFast(DEBUG_PIN_3, debug_state_3=!debug_state_3);
#else
#define DEBUG_SET_3
#define DEBUG_CLEAR_3
#define DEBUG_TOG_3
#endif

#ifdef DEBUG_PIN_4
extern int debug_state_4;
#define DEBUG_SET_4 digitalWriteFast(DEBUG_PIN_4, HIGH);
#define DEBUG_CLEAR_4 digitalWriteFast(DEBUG_PIN_4, LOW);
#define DEBUG_TOG_4 digitalWriteFast(DEBUG_PIN_4, debug_state_4=!debug_state_4);
#else
#define DEBUG_SET_4
#define DEBUG_CLEAR_4
#define DEBUG_TOG_4
#endif

#endif //TRJ_ENV_HOST

void ser_printf(const char* fmt, ...);