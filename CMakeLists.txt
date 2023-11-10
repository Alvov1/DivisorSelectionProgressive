cmake_minimum_required(VERSION 3.22)
project(DivisorSelectionProgressive CUDA CXX)

set(CMAKE_CUDA_STANDARD 23)

include(FetchContent)
FetchContent_Declare(AesiMultiprecision
        GIT_REPOSITORY https://github.com/Alvov1/Aesi-Multiprecision.git
        GIT_TAG main
)
FetchContent_MakeAvailable(AesiMultiprecision)

add_executable(DivisorSelectionProgressive main.cu)
target_include_directories(DivisorSelectionProgressive PRIVATE ${AesiMultiprecision_SOURCE_DIR})

