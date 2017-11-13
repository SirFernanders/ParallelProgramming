cmake_minimum_required(VERSION 3.8)
project(ParallelProgramming)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        pth_pi_mutex ipp-source-use/ch4/pth_pi_mutex.c)

add_executable(pth_pi_mutex ipp-source-use/ch4/pth_pi_mutex.c)
set_target_properties(pth_pi_mutex PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/ipp-source-use/ch4/")