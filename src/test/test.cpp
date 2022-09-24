#include <iostream>
#include <sstream>
#include <string>
#include <catch2/catch_test_macros.hpp>
#define CATCH_CONFIG_MAIN

using namespace std;

std::vector<int> extractIntegerWords(const string& str)
{
    stringstream ss(str);
    std::vector<int> result;

    string temp;
    int found;
    while (!ss.eof()) {

        ss >> temp;

        if (stringstream(temp) >> found)
            result.push_back(found);

        temp = "";
    }

    return result;
}