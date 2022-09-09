#pragma once

#include <Arduino.h>
#include <limits.h>
#include "trj_const.h"

// Configuration record for one axis
// 8 Bytes
struct AxisConfig {

    uint8_t axis;           // Axis number

    uint8_t step_pin;       // Step output, or quadture b
    uint8_t direction_pin;  // Direction output, or quadrature b
    uint8_t enable_pin;

    uint8_t step_high_value; // Whether step is HIGH or LOW when enabled. 
    uint8_t direction_high_value; 
    uint8_t enable_high_value; 

    uint8_t step_output_mode; // OUTPUT or OUTPUT_OPEN_DRAIN
    uint8_t direction_output_mode;
    uint8_t enable_output_mode;

    uint8_t pad1;
    uint8_t pad2;

    uint32_t v_max;
    uint32_t a_max;
};

// Main Configuration class
struct Config {

    uint8_t n_axes = 0;         // Number of axes
    uint8_t interrupt_delay=INTERRUPT_DELAY;    // How often interrupt is called, in microseconds
    uint8_t segment_complete_pin=0; // Pin on which to signal that a segment is complete
    uint8_t limit_pin=0; // Pin to recieve signals that the encoder foind a limit
    bool debug_print = true;
    bool debug_tick = true;
};
