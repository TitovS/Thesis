cmake_minimum_required(VERSION 3.2)
project(Thesis)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(Thesis ${SOURCE_FILES} Stat.cpp Stat.h Lorenz.cpp Lorenz.h 1.txt)