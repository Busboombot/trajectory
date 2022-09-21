cmake_minimum_required(VERSION 3.23)
project(pyplanner LANGUAGES CXX VERSION 0.1)
set (CMAKE_CXX_STANDARD 11)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(PythonInterp 3)
find_package(PythonLibs 3 REQUIRED)
find_package(Boost COMPONENTS system python3 REQUIRED)


python3_add_library(pyplan  test.cpp)



target_link_libraries(pyplan PUBLIC
        ${Boost_LIBRARIES}
        libplanner
        ${PYTHON_LIBRARIES}

        )
target_include_directories(pyplan PUBLIC
        ${Boost_INCLUDE_DIRS}
        ${PYTHON_INCLUDE_DIRS}
        )

#set_target_properties(pyplan PROPERTIES VERSION ${PROJECT_VERSION})
#set_target_properties(pyplan PROPERTIES SOVERSION ${PROJECT_VERSION})

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")