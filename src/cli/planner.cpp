

#include <string>
#include <iostream>
#include <vector>
#include <sstream>

#include "trj_joint.h"
#include "trj_planner.h"

#include <chrono>
using namespace std;

std::vector<int> extractIntegerWords(string str)
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

using Ints = vector<int>;
using Moves = vector<Ints>;

int main() {

    vector<Joint> joints;
    Moves moves;
    int n_joints;
    int line_n = 0;
    Planner *planner = nullptr;

    for (std::string line; std::getline(std::cin, line);) {
        //cout << line_n << " " << line << endl;
        Ints ints = extractIntegerWords(line);
        if (line_n == 0) {
            n_joints = ints[0];
        } else if (line_n <= n_joints) {
            joints.emplace_back(line_n-1, ints[0],ints[1]);
        } else {
            moves.push_back(ints);
        }
        line_n += 1;
    }

    auto start = chrono::steady_clock::now();
    planner = new Planner(joints);
    for(Ints &m : moves) {
        planner->move(m);
    }
    auto end = chrono::steady_clock::now();
    auto diff = chrono::duration_cast<chrono::nanoseconds>(end - start);

    json j = planner->dump();
    j["_time"] = diff.count();

    cout << j << endl;

    return 0;
}

