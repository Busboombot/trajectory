#include <boost/python.hpp>

#include "trj_block.h"
#include "trj_segment.h"

using namespace boost::python;

char const* greet()
{
    return "hello, world";
}


BOOST_PYTHON_MODULE(pyplan)
{
    using namespace boost::python;
    def("greet", greet);


}


