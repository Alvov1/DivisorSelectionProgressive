cmake_minimum_required(VERSION 3.25)
project(PrimeGenerator CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

find_package(primesieve REQUIRED)

add_executable(PrimeGenerator main.cpp)
target_link_libraries(PrimeGenerator primesieve::primesieve)

