cmake_minimum_required(VERSION 3.5)
project(solver_clion)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")

# find_package(Eigen3 3.3 REQUIRED)
#include_directories(${EIGEN3_INCLUDE_DIR})

set(EIGEN_DIR "/usr/local/include/eigen3")
include_directories(${EIGEN_DIR})

add_executable(solver_clion
        main.cpp
        equations.cpp
        equations.hpp
        parity.cpp
        parity.hpp
        constants.hpp
        eigenvalue.cpp
        eigenvalue.hpp
        precise_eigenvalue_finder.cpp
        precise_eigenvalue_finder.hpp
        preliminary_eigenvalue_finder.cpp
        preliminary_eigenvalue_finder.hpp
        )

#target_link_libraries(solver_clion Eigen3::Eigen)
target_link_libraries(solver_clion ${EIGEN_DIR})

find_package(GSL REQUIRED)
target_link_libraries(solver_clion GSL::gsl GSL::gslcblas)
