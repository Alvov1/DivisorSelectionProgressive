#include <iostream>
#include <fstream>
#include <filesystem>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/device_ptr.h>
#include "AesiMultiprecision.h"

std::vector<uint64_t> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::exists(fromLocation) || !std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");

    std::ifstream input(fromLocation, std::ios::binary);
    uint64_t buffer {}; input.read(reinterpret_cast<char*>(&buffer), sizeof(uint64_t));

    std::vector<uint64_t> primes (buffer);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(uint64_t));

    return primes;
}

//gpu void kernel(const Aesi& n, const thrust::device_vector<uint64_t>& primes) {
//    const auto threadIdx = blockDim.x * blockIdx.x + threadIdx.x,
//        threads = gridDim.x * blockDim.x,
//        max_it = 400000 / threads,
//        bStart = 2 + blockIdx.x,
//        bInc = gridDim.x,
//        B_MAX = 2000000000;
//
//    const auto checkWriteRepeat = [&n] (const Aesi<512>& value) {
//        if(value < 2 || value >= n) return false;
//        char buffer[512] {}; value.getString<10>(buffer, 512);
//        return printf("Found divisor: %s\n", buffer) > 14;
//    };
//
//    Aesi a = threadIdx * max_it + 2, e = 1;
//    for(unsigned B = bStart; B < B_MAX; B += bInc) {
//        const auto primeUl = primes[0];
//
//        for(unsigned pi = 0; primeUl < B; ++pi) {
//            const unsigned power = log(static_cast<double>(B)) / log(static_cast<double>(primeUl));
//            e *= static_cast<uint64_t>(pow(static_cast<double>(primeUl), static_cast<double>(power)));
//            primeUl = primes[pi + 1];
//        }
//
//        if(e == 1) continue;
//
//        for(unsigned it = 0; it < max_it; ++it) {
//            if(checkWriteRepeat(Aesi<512>::gcd(a, n)))
//                return;
//
//            if(checkWriteRepeat(Aesi<512>::gcd(Aesi<512>::powm(a, e, n) - 1, n)))
//                return;
//
//            a += threads * max_it;
//        }
//    }
//}

__global__ void kernel(Aesi<512>* const numbers) {
    const auto tid = blockDim.x * blockIdx.x + threadIdx.x;
    if(tid > 0) return;

    Aesi left = *numbers;
    numbers[1] = left * 10 + left ^ "12345678901234567890123456789012345678901234567890";
}

int main(int argc, const char* const* const argv) {
    if(argc < 4)
        return std::printf("Usage: %s factorize <number> <primes location>", argv[0]);

    Aesi<512> number = std::string_view(argv[2]);
    std::cout << "Factorizing number " << std::hex << std::showbase << number << '.' << std::endl;

//    const thrust::device_vector<uint64_t> primeTable = loadPrimes(argv[3]);
//    std::cout << "Loaded prime table of " << primeTable.size() << " elements." << std::endl;

    thrust::device_vector<Aesi<512>> numbers = { number, {} };

    kernel<<<32, 32>>>(thrust::raw_pointer_cast(numbers.data())); const auto code = cudaDeviceSynchronize();
    if (code != cudaSuccess)
        return std::printf("Kernel launch failed: %s.\n", cudaGetErrorString(code));

    std::cout << "Kernel completed. Got product: " << numbers[1] << std::endl;

    return 0;
}
