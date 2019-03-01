cmake_minimum_required(VERSION 3.13)
project(normalize)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED  11)

set(EIGEN_DIR "/usr/local/include/eigen3")
include_directories(${EIGEN_DIR})


add_executable(normalize main.cpp
                         channel.cpp
                         channel.hpp
                         wavefunction.cpp
                         wavefunction.hpp
                         filelister.cpp
                         filelister.hpp)

target_link_libraries(normalize ${EIGEN_DIR})

find_package(GSL REQUIRED)
target_link_libraries(normalize GSL::gslcblas GSL::gsl)
