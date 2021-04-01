cmake_minimum_required(VERSION 3.16)
project(ss)
include_directories(include)
set(CMAKE_CXX_STANDARD 17)

add_executable(ss src/main.cpp src/schreier.cpp include/schreier.h)