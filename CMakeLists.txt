cmake_minimum_required(VERSION 3.12)
project(MPM)

set(CMAKE_CXX_STANDARD 14)
find_package(Eigen3 REQUIRED)
add_executable(MPM main.cpp MPM.h)
target_link_libraries (MPM Eigen3::Eigen)
