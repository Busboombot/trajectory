cmake_minimum_required(VERSION 3.24)
set (CMAKE_CXX_STANDARD 11)

project (trajectory VERSION 0.0.1)
 
add_subdirectory(planner) 


include(CTest)
find_package(Catch2 3 REQUIRED)
include(Catch)

# We check if this is the main file
# you don't usually want users of your library to
# execute tests as part of their build
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    add_subdirectory(test)
endif ()

enable_testing()
