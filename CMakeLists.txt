cmake_minimum_required(VERSION 3.10.0 FATAL_ERROR)
project(PARTIO_EXAMPLE)
set (CMAKE_CXX_STANDARD 11)
find_package(Eigen3 REQUIRED)
set(EXTERNAL ${CMAKE_BINARY_DIR}/external)

include(partio.cmake)

file(GLOB SRC "*.h" "*.cpp")

#add_executable(main test.cpp)
#target_link_libraries(main partio)
#add_dependencies(main partio-dl)

#add_executable(main main.cpp Particle.cpp Particle.h Grid.h Grid.cpp Constants.h)
#target_link_libraries(main partio Eigen3::Eigen)
add_executable(main test.cpp)
target_link_libraries(main partio)
add_dependencies(main partio-dl)