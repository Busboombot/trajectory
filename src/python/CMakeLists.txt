
cmake_minimum_required(VERSION 3.23)

project(pyplanner LANGUAGES CXX VERSION 0.1)
set (CMAKE_CXX_STANDARD 11)

set(CMAKE_SHARED_MODULE_PREFIX "")

add_library(pyplanlib SHARED
    test.cpp
        )

PYTHON_ADD_MODULE(pyplan test.cpp)

find_package(Boost COMPONENTS system python3 REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

target_link_libraries(pyplanlib PUBLIC
        ${Boost_LIBRARIES}
        planner
        ${PYTHON_LIBRARIES}

        )
target_include_directories(pyplanlib PUBLIC
        ${Boost_INCLUDE_DIRS}
        ${PYTHON_INCLUDE_DIRS}
      
        )


set_target_properties(pyplanlib PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(pyplanlib PROPERTIES SOVERSION ${PROJECT_VERSION})

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")