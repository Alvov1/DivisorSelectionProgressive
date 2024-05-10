cmake_minimum_required(VERSION 3.25)
project(DivisorSelectionProgressive CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(AesiMultiprecision
        GIT_REPOSITORY https://github.com/Alvov1/Aesi-Multiprecision.git
        GIT_TAG split-signed-unsigned-integers
)
FetchContent_MakeAvailable(AesiMultiprecision)

add_executable(PowersGeneration generate.cpp)
target_include_directories(PowersGeneration PRIVATE ${AesiMultiprecision_SOURCE_DIR})
