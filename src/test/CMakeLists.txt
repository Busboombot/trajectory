
cmake_minimum_required(VERSION 3.23)

project(planner_tests LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 11)

find_package(Boost COMPONENTS system filesystem REQUIRED)

set( SOURCES_TESTS
        test.cpp
        test_basic.cpp
        test_block.cpp
        test_segment.cpp
        test_planner.cpp
        test_stepper.cpp

)

add_executable(test_planner ${SOURCES_TESTS} )

target_link_libraries(test_planner PRIVATE Catch2::Catch2WithMain libplanner ${Boost_LIBRARIES})

target_include_directories(test_planner PUBLIC
        ${Boost_INCLUDE_DIRS}
        )

catch_discover_tests(test_planner)
