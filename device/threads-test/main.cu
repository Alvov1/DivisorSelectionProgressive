#include <iostream>
#include <fstream>
#include <filesystem>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

__global__ void kernel(unsigned* primes, std::size_t primesCount) {
    const unsigned threadId = blockDim.x * blockIdx.x + threadIdx.x,
      threadsCount = gridDim.x * blockDim.x;
    if(threadId % 32768 != 0) return;

    if(threadId != 0)
      printf("Thread %u. Prime: %u\n", threadId, primes[threadId % primesCount]);
    else
      printf("Thread %u. Threads count %u. Prime: %u\n", threadId, threadsCount, primes[threadId % primesCount]);
}

using primeType = unsigned;
thrust::host_vector<primeType> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");

    const auto primesCount = std::filesystem::file_size(fromLocation) / sizeof(primeType);
    std::ifstream input(fromLocation, std::ios::binary);

    std::vector<primeType> primes (primesCount);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(primeType));

    return primes;
}

int main(int argc, const char* const* const argv) {
    try {
        if (argc < 4)
            return std::printf("Usage: %s <primes location> <lThreads> <rThreads>.\n", argv[0]);

        const std::filesystem::path primesLocation (argv[1]);
        thrust::device_vector<primeType> primes = loadPrimes(primesLocation);
        std::printf("Loaded prime table with %lu elements\n", primes.size());

        const unsigned lThreads = std::stoi(argv[2]), rThreads = std::stoi(argv[3]);
        std::printf("Launching kernel <<<%u, %u>>>.\n", lThreads, rThreads);

        kernel<<<lThreads, rThreads>>>(
          thrust::raw_pointer_cast(primes.data()),
          primes.size());

        const auto code = cudaDeviceSynchronize();
        if (code != cudaSuccess)
            return std::printf("Kernel launch failed: %s.\n", cudaGetErrorString(code));

        std::printf("Kernel completed.\n");
    } catch (const std::exception& e) {
        return std::printf("Failed: %s.\n", e.what());
    }
}