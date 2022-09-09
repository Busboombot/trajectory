cmake_minimum_required(VERSION 3.24)
set (CMAKE_CXX_STANDARD 11)

project (trajectory VERSION 0.0.1)
 
add_subdirectory(planner) 
include(CTest)


# We check if this is the main file
# you don't usually want users of your library to
# execute tests as part of their build
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    add_subdirectory(test)
endif ()

enable_testing()

#add_executable(hello-cmake hello.cpp)

target_include_directories(planner PUBLIC planner )
target_compile_definitions(planner PUBLIC "TRJ_ENV_HOST")
target_compile_definitions(planner PUBLIC INTERRUPT_DELAY=4)

