
cmake_minimum_required(VERSION 3.24)

project(planner)

add_library( ${PROJECT_NAME}
trj_jointss.cpp
trj_move.cpp
trj_planner.cpp
trj_segment.cpp
trj_stepdriver.cpp
trj_util.cpp
)

target_include_directories( ${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}
)