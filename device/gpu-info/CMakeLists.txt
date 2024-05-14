cmake_minimum_required(VERSION 3.22)
project(DivisorSelectionProgressive CUDA CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(CudaDeviceInfo
        GIT_REPOSITORY https://github.com/Alvov1/CudaDeviceInfo.git
        GIT_TAG main
)
FetchContent_MakeAvailable(CudaDeviceInfo)

add_executable(GetGpuInfo main.cu)
set_property(TARGET GetGpuInfo PROPERTY CUDA_ARCHITECTURES 75)
target_include_directories(GetGpuInfo PRIVATE ${CudaDeviceInfo_SOURCE_DIR})