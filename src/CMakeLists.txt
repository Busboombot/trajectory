cmake_minimum_required(VERSION 3.23)
set (CMAKE_CXX_STANDARD 11)

project (trajectory VERSION 0.0.1)


include(CTest)
find_package(Catch2 3 REQUIRED)
include(Catch)


find_package(Boost 1.79.0 REQUIRED COMPONENTS
        filesystem system python310 )

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find boost!")
endif()

# We check if this is the main file
# you don't usually want users of your library to
# execute tests as part of their build
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    add_subdirectory(test)
endif ()

add_subdirectory(planner)
add_subdirectory(cli)

enable_testing()
