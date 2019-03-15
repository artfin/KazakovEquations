cmake_minimum_required(VERSION 3.5)
project(offdiag_coriolis)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED  11)

set(EIGEN_DIR "/usr/local/include/eigen3")
include_directories(${EIGEN_DIR})

add_executable(offdiag_coriolis
               main.cpp
               channel.cpp
               channel.hpp
               wavefunction.cpp
               wavefunction.hpp
               filelister.cpp
               filelister.hpp)

find_package(GSL REQUIRED)
target_link_libraries(offdiag_coriolis GSL::gslcblas GSL::gsl)
