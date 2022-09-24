
cmake_minimum_required(VERSION 3.23)

project(cli LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 11)

find_package(Boost COMPONENTS program_options REQUIRED)

add_executable(planner planner.cpp )

target_link_libraries(planner PUBLIC
        ${Boost_LIBRARIES}
        libplanner
        )
target_include_directories(planner PUBLIC
        ${Boost_INCLUDE_DIRS}
        )

add_executable(stepper stepper.cpp )
target_link_libraries(stepper PRIVATE  libplanner)