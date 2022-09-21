
cmake_minimum_required(VERSION 3.23)

project(cli LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 11)


add_executable(planner planner.cpp )
target_link_libraries(planner PRIVATE  libplanner)

add_executable(stepper stepper.cpp )
target_link_libraries(stepper PRIVATE  libplanner)