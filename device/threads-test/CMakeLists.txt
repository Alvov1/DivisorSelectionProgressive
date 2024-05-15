cmake_minimum_required(VERSION 3.25)
project(DivisorSelectionProgressive CUDA CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

add_executable(PrimesSelection main.cu)
#set_property(TARGET PrimesSelection PROPERTY CUDA_ARCHITECTURES 75)