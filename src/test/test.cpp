#include <iostream>

#include <catch2/catch_test_macros.hpp>
#define CATCH_CONFIG_MAIN

TEST_CASE("Is this thing on?", "[basic]")
{
    std::cout << "This is test.cpp" << std::endl;
   
    REQUIRE(true);
}