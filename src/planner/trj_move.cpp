

#include <cstdint> 
#include <iostream>
#include <iomanip>

#include "trj_move.h"

using std::array;
using std::cout;
using std::endl;
using std::setw;
using std::left;
using std::right;
using std::ostream;

ostream &operator<<( ostream &output, const Move &m ) {
    
    output << "[Move #" << m.seq << " " << (int)m.move_type <<  " t=" <<m.t << " (" ;

    for(auto &xi : m.x){
        output << xi << ", ";
    }

    output << ")]" << endl; 

    return output;
}
