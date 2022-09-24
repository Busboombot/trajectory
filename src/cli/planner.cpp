

#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <exception>

#include <boost/program_options.hpp>

#include "trj_joint.h"
#include "trj_planner.h"

#include <chrono>
using namespace std;
namespace po = boost::program_options;

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

void loadData(vector<Joint> &joints, Moves &moves){

    int n_joints;
    int line_n = 0;


    for (std::string line; std::getline(std::cin, line);) {
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

}

Planner *makePlanner(vector<Joint> &joints, Moves &moves){

    Planner *planner = new Planner(joints);
    for(Ints &m : moves) {
        planner->move(m);
    }

    return planner;
}

int main(int ac, char **av) {

    vector<Joint> joints;
    Moves moves;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce help message")
            ("planner,p",  "Load moves into the planner and print it. ")
            ("stepper,s",  "Load moves into the planner run steppers ")
            ("json,j",  "Output JSON ")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << endl;
        return 1;
    }

    loadData(joints, moves);

    auto start = chrono::steady_clock::now();
    Planner *planner = makePlanner(joints, moves);
    auto end = chrono::steady_clock::now();
    auto diff = chrono::duration_cast<chrono::microseconds>(end - start);

    if (!vm.count("json")){
        cout << "Processed moves in " << diff.count() << "μs" << endl;
    }

    if(vm.count("planner")){
        if (vm.count("json")){
            json j = planner->dump();
            j["_time"] = diff.count();
            cout << j << endl;
        } else {
            cout << *planner << endl;
        }
    } else if(vm.count("stepper")){
        if (vm.count("json")){

        } else {
            cout << *planner << endl;
        }
    }


    return 0;
}

