
cmake_minimum_required(VERSION 3.23)

project(planner_tests LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 11)


set( SOURCES_TESTS
        test.cpp
        test_basic.cpp
        test_block.cpp
        test_planner.cpp
)

add_executable(test_planner ${SOURCES_TESTS} )

target_link_libraries(test_planner PRIVATE Catch2::Catch2WithMain planner)

catch_discover_tests(test_planner)

