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
FetchContent_Declare(STDStreamTimer
        GIT_REPOSITORY https://github.com/Alvov1/STD-Stream-Timer.git
        GIT_TAG main
)
FetchContent_MakeAvailable(STDStreamTimer)

add_executable(DivisorPrimes primes.cpp)
target_include_directories(DivisorPrimes PRIVATE ${AesiMultiprecision_SOURCE_DIR} ${STDStreamTimer_SOURCE_DIR})

add_executable(DivisorPowers powers.cpp
        ../tmain.cu)
target_include_directories(DivisorPowers PRIVATE ${AesiMultiprecision_SOURCE_DIR} ${STDStreamTimer_SOURCE_DIR})

add_executable(DivisorPreparations ../powers/generate.cpp)
target_include_directories(DivisorPreparations PRIVATE ${AesiMultiprecision_SOURCE_DIR} ${STDStreamTimer_SOURCE_DIR})

add_executable(CircularBuffer main.cpp)
target_include_directories(CircularBuffer PRIVATE ${AesiMultiprecision_SOURCE_DIR} ${STDStreamTimer_SOURCE_DIR})