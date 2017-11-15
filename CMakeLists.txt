cmake_minimum_required(VERSION 3.8)
project(ParallelProgramming)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        pth_pi_mutex ipp-source-use/ch4/pth_pi_mutex.c)

add_executable(pth_pi_mutex ipp-source-use/ch4/pth_pi_mutex.c)
set_target_properties(pth_pi_mutex PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/ipp-source-use/ch4/")
add_executable(4.1_pth_pi_mutex "Assignment 4/4.1_pth_pi_mutex.c")
set_target_properties(4.1_pth_pi_mutex PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/Assignment 4/")
add_executable(4.2_pth_pi_mutex "Assignment 4/4.2_pth_pi_mutex.c")
set_target_properties(4.2_pth_pi_mutex PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/Assignment 4/")
add_executable(4_3 "Assignment 4/4_3.c")
set_target_properties(4_3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/Assignment 4/")
add_executable(Sp "Assignment 4/Sp.c")
set_target_properties(Sp PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/Assignment 4/")