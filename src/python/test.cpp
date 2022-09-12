#include <boost/python.hpp>

#include "trj_block.h"

char const* greet()
{
    return "hello, world";
}


BOOST_PYTHON_MODULE(pyplan)
{
    using namespace boost::python;
    def("greet", greet);
}
