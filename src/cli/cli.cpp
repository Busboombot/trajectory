

#include <string>
#include <iostream>
#include <vector>
#include <sstream>

#include "trj_joint.h"
#include "trj_planner.h"

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


int main() {

    vector<Joint> joints;
    int n_joints;
    int line_n = 0;
    Planner *planner = nullptr;

    for (std::string line; std::getline(std::cin, line);) {
        //cout << line_n << " " << line << endl;
        auto ints = extractIntegerWords(line);
        if (line_n == 0) {
            n_joints = ints[0];
        } else if (line_n <= n_joints) {
            joints.emplace_back(line_n-1, ints[0],ints[1]);
        } else {

            if (planner == nullptr){
                planner = new Planner(joints);
            }
            planner->move(ints);
        }
        line_n += 1;
    }

    cout << planner->dump() << endl;

    return 0;
}

