cmake_minimum_required(VERSION 3.22)
project(DivisorSelectionProgressive CUDA)

include(FetchContent)
FetchContent_Declare(AesiMultiprecision
        GIT_REPOSITORY https://github.com/Alvov1/Aesi-Multiprecision.git
        GIT_TAG main
)
FetchContent_MakeAvailable(AesiMultiprecision)

add_executable(DivisorSelectionProgressive main.cpp)
set_property(TARGET DivisorSelectionProgressive PROPERTY CUDA_STANDARD 20)
target_include_directories(DivisorSelectionProgressive PRIVATE ${AesiMultiprecision_SOURCE_DIR})

