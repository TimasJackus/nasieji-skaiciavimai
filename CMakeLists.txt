cmake_minimum_required(VERSION 3.15)
project(nasieji)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-fopenmp")
add_executable(nasieji main.cpp)