cmake_minimum_required(VERSION 3.9)

project(kludemo)

include_directories(/home/liam/repos/SuiteSparse-5.10.1/include)

link_directories(/home/liam/repos/SuiteSparse-5.10.1/lib)

add_definitions(-O3 -Wall -std=c99)

add_executable(kludemo kludemo.c)

target_link_libraries(kludemo cholmod amd colamd camd ccolamd metis lapack blas klu m)
