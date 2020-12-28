cmake_minimum_required(VERSION 3.13.4)
project(nasieji)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)
add_executable(nasieji main.cpp)
target_link_libraries(nasieji PRIVATE OpenMP::OpenMP_CXX)
