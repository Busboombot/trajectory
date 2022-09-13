// Copyright (c) 2018 Daniel Lovasko.
// All Rights Reserved
//
// Distributed under the terms of the 2-clause BSD License. The full
// license is in the file LICENSE, distributed as part of this software.
#pragma once
#include <string>
#include <sstream>
#include <ostream>

/// ANSI coloring for strings.
namespace col {

  /// Color hues.
  enum hue {
    black, red, green, yellow, hblue, blue, magenta, cyan, light_gray, dark_gray,
    light_red, light_green, light_yellow, light_blue, light_magenta,
    light_cyan, white, def
  };

  /// Retrieve the ANSI color code if used in context of a foreground.
  int fg_code(const hue h) {
    if (h == black)         return 30;
    if (h == red)           return 31;
    if (h == green)         return 32;
    if (h == yellow)        return 33;
    if (h == blue)          return 34;
    if (h == magenta)       return 35;
    if (h == cyan)          return 36;
    if (h == light_gray)    return 37;
    if (h == dark_gray)     return 90;
    if (h == light_red)     return 91;
    if (h == light_green)   return 92;
    if (h == light_yellow)  return 93;
    if (h == light_blue)    return 94;
    if (h == light_magenta) return 95;
    if (h == light_cyan)    return 96;
    if (h == white)         return 97;

    return 39; // def color
  }

  /// Retrieve the ANSI color code if used in context of a background.
  int bg_code(const hue h) {
    if (h == black)         return 40;
    if (h == red)           return 41;
    if (h == green)         return 42;
    if (h == yellow)        return 43;
    if (h == blue)          return 44;
    if (h == hblue)         return 104;
    if (h == magenta)       return 45;
    if (h == cyan)          return 46;
    if (h == light_gray)    return 47;
    if (h == dark_gray)     return 100;
    if (h == light_red)     return 101;
    if (h == light_green)   return 102;
    if (h == light_yellow)  return 103;
    if (h == light_blue)    return 104;
    if (h == light_magenta) return 105;
    if (h == light_cyan)    return 106;
    if (h == white)         return 107;

    return 49; // def color
  }

  /// Make a color based on selected attributes.
  std::string make(const hue fg, const hue bg, const bool bold,
                   const bool uline, const bool blink) {
    std::stringstream ss;

    ss << "\x1b[";
    ss << fg_code(fg) << ";";
    ss << bg_code(bg);
    if (bold)  ss << ";" << 1;
    if (uline) ss << ";" << 4;
    if (blink) ss << ";" << 5;
    ss << "m";

    return ss.str();
  }

  /// Output stream manipulator for resetting the color mode.
  std::ostream& reset(std::ostream& os) {
    os << "\x1b[0m";
    return os;
  }
}
