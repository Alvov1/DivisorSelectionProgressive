cmake_minimum_required(VERSION 3.25)
project(DivisorSelectionProgressive CUDA CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)

include(FetchContent)
FetchContent_Declare(AesiMultiprecision
        GIT_REPOSITORY https://github.com/Alvov1/Aesi-Multiprecision.git
        GIT_TAG main
)
FetchContent_MakeAvailable(AesiMultiprecision)

add_executable(DivisorSelectionProgressive main.cu)
target_include_directories(DivisorSelectionProgressive PRIVATE ${AesiMultiprecision_SOURCE_DIR})

