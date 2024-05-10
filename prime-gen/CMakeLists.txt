cmake_minimum_required(VERSION 3.22)
project(PrimeGenerator CXX)

set(CMAKE_CXX_STANDARD 20)
find_package(primesieve REQUIRED)

add_executable(PrimeGenerator main.cpp)
target_link_libraries(PrimeGenerator primesieve::primesieve)
target_include_directories(PrimeGenerator PRIVATE ${STDStreamTimer_SOURCE_DIR})
