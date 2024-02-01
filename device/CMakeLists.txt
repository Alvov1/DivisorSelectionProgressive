cmake_minimum_required(VERSION 3.25)
project(DivisorSelectionProgressive CUDA CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)


include(FetchContent)
FetchContent_Declare(AesiMultiprecision
        GIT_REPOSITORY https://github.com/Alvov1/Aesi-Multiprecision.git
        GIT_TAG main
)
FetchContent_MakeAvailable(AesiMultiprecision)
FetchContent_Declare(STDStreamTimer
        GIT_REPOSITORY https://github.com/Alvov1/STD-Stream-Timer.git
        GIT_TAG main
)
FetchContent_MakeAvailable(STDStreamTimer)

add_executable(DivisorPrimes primes.cu)
set_property(TARGET DivisorPrimes PROPERTY CUDA_ARCHITECTURES 75)
target_include_directories(DivisorPrimes PRIVATE ${AesiMultiprecision_SOURCE_DIR} ${STDStreamTimer_SOURCE_DIR})

add_executable(DivisorPowers powers.cu)
set_property(TARGET DivisorPowers PROPERTY CUDA_ARCHITECTURES 75)
target_include_directories(DivisorPowers PRIVATE ${AesiMultiprecision_SOURCE_DIR} ${STDStreamTimer_SOURCE_DIR})

