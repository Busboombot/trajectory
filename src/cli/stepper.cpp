/* Read integers in text from stdin, construct steppers
 * and step out steps.
 *
 * Each row of the input is one stepper phase for multiple axes, with
 * the line consisting of 1 or more groups of three numbers, eah representing
 * (x, vi, vf). Each group of three rows represents the 3 phases for a
 * blocks, A, C, D. So, each row where line_no%3 == 0 is an A phase,
 * line_no%3 == 1 is a C phase, and line_no%3 == 2 is a D phase.
 *
 * */
#include <iostream>
#include <vector>
#include <sstream>
#include <array>
#include <list>

#include "trj_joint.h"
#include "trj_planner.h"
#include "trj_stepper.h"

#include <chrono>
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

using Ints = vector<int>; // One input line of integers. Ought to be divisible by  3

using StepperBlock = array<StepperPhase, 3>;
using Blocks = vector<StepperBlock>;
using Segments = vector<Blocks>;

ostream &operator<<(ostream &output, const StepperPhase &s) {
    output << "(" << s.vi << "," << s.x << "," << s.vf << ")";
    return output;
}

ostream &operator<<(ostream &output, const StepperBlock &sb) {
    output << "[" << sb[0] << "/" << sb[1] << "/" << sb[2] << "]";
    return output;
}

ostream &operator<<(ostream &output, const Blocks &b) {
   for(const StepperBlock &sb: b){
       output << sb << " ";
   }
   return output;
}

ostream &operator<<(ostream &output, const Segments &s) {
    for(const Blocks &b: s) {
        output << b << endl;
    }
    return output;
}


int main() {

    Segments segments;
    double dtime = 5./1e6; // 5 us

    /// Load all of the lines in to vectors

    int line_n = 0;
    int n_axes = 0;
    for (std::string line; std::getline(std::cin, line);) {
        if (line[0] == ' ' || line[0] == '#') continue;
        Ints ints = extractIntegerWords(line); // Get all integers on a line
        if (ints.size() % 3 != 0) continue; // Each group of 3 is a block

        if (line_n%3 == 0) {
            n_axes = ints.size()/3;
            segments.push_back(Blocks(n_axes));
        }

        for(int i = 0; i < ints.size(); i+= 3){
            segments.back()[i/3][line_n%3] = StepperPhase{ints[i],(float)ints[i+1], (float)ints[i+2]};
        }

        line_n += 1;
    }

    cout << segments << endl;

    auto steppers = vector<StepperState>(n_axes);

    int seg_n = 0;
    for (Blocks &b: segments) {

        for(int i =0; i < steppers.size(); i++){
            steppers[i].loadPhases(b[i]);
        }

        int doneCount = 0;
        int step_n = 0;
        auto dist = vector<int>(n_axes);

        do {
            int sn = 0;
            for( auto &ss: steppers){
                doneCount += ss.isDone();
                dist[sn++] += ss.next(dtime);
            }
            step_n += 1;
        } while(doneCount == 0);

        cout << seg_n << " " <<  step_n << " " ;
        for(auto i: dist) cout << i << " ";
        cout << endl;

        seg_n++;
    }



    return 0;
}
