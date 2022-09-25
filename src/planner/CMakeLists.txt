
cmake_minimum_required(VERSION 3.23)
set (CMAKE_CXX_STANDARD 11)

project(libplanner VERSION 0.1)

add_library( ${PROJECT_NAME}
        trj_move.cpp
        trj_segment.cpp
        trj_util.cpp
        trj_block.cpp
        trj_joint.cpp
        trj_planner.cpp
        trj_stepper.cpp trj_stepper.h trj_segstepper.cpp)

configure_file(version.h.in version.h)
target_compile_definitions(${PROJECT_NAME} PUBLIC "TRJ_ENV_HOST")
target_compile_definitions(${PROJECT_NAME} PUBLIC INTERRUPT_DELAY=4)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR} )