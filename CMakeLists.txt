cmake_minimum_required(VERSION 3.12)
project(MPM)

set(CMAKE_CXX_STANDARD 14)
find_package(Eigen3 REQUIRED)
include_directories(/usr/local/include)
add_executable(MPM main.cpp MPM.h MPM.cpp)
target_link_libraries (MPM Eigen3::Eigen)
