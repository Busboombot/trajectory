
cmake_minimum_required(VERSION 3.23)

project(trj_planner LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 11)

set( SOURCES_TESTS
        cli.cpp)

add_executable(trj_planner ${SOURCES_TESTS} )

target_link_libraries(trj_planner PRIVATE  planner)

