cmake_minimum_required(VERSION 3.24)

project (trajectory VERSION 0.0.1)
 
add_subdirectory(planner) 

#add_executable(hello-cmake hello.cpp)

target_include_directories(trajectory PUBLIC planner )
