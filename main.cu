#include <iostream>
#include <fstream>
#include <filesystem>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/device_ptr.h>
#include "AesiMultiprecision.h"
#include "Timer.h"

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

__global__ void kernel(const Aesi<512>* numberAndFactor, const uint64_t* const primes, std::size_t primesCount) {
    const auto threadId = blockDim.x * blockIdx.x + threadIdx.x,
            threads = gridDim.x * blockDim.x,
            max_it = 400000 / threads,
            bStart = 2 + blockIdx.x,
            bInc = gridDim.x,
            B_MAX = 2000000000U;

    const Aesi n = numberAndFactor[0], factor = numberAndFactor[1];
    const auto checkWriteRepeat = [&n](const Aesi<512> &value) {
        if (value < 2 || value >= n) return false;
        char buffer[512]{};
        value.getString<10>(buffer, 512);
        return printf("Found divisor: %s\n", buffer) > 14;
    };
    if(!factor.isZero()) return;

    Aesi a = threadId * max_it + 2, e = 1;
    for (unsigned B = bStart; B < B_MAX; B += bInc) {
        auto primeUl = primes[0];

        for (unsigned pi = 0; primeUl < B; ++pi) {
            if(!factor.isZero()) return;
            const unsigned power = log(static_cast<double>(B)) / log(static_cast<double>(primeUl));
            e *= static_cast<uint64_t>(pow(static_cast<double>(primeUl), static_cast<double>(power)));
            primeUl = primes[pi + 1];
        }

        if (e == 1) continue;

        for (unsigned it = 0; it < max_it; ++it) {
            if(!factor.isZero()) return;
            if (checkWriteRepeat(Aesi<512>::gcd(a, n)))
                return;

            if (checkWriteRepeat(Aesi<512>::gcd(Aesi<512>::powm(a, e, n) - 1, n)))
                return;

            a += threads * max_it;
        }
    }
}

int main(int argc, const char* const* const argv) {
    if(argc < 4)
        return std::printf("Usage: %s factorize <number> <primes location>", argv[0]);

    thrust::device_vector<Aesi<512>> numberAndFactor = { { std::string_view(argv[2]) }, {} };
    Timer::init() << "Factorizing number " << std::hex << std::showbase << numberAndFactor[0] << '.' << Timer::endl;

    const thrust::device_vector<uint64_t> primes = loadPrimes(argv[3]);
    Timer::out << "Loaded prime table of " << primes.size() << " elements." << Timer::endl;

    kernel<<<32, 32>>>(
            thrust::raw_pointer_cast(numberAndFactor.data()),
            thrust::raw_pointer_cast(primes.data()),
            primes.size());

    const auto code = cudaDeviceSynchronize();
    if (code != cudaSuccess)
        return std::printf("Kernel launch failed: %s.\n", cudaGetErrorString(code));
    Timer::out << "Kernel completed. Founded factor: " << std::hex << std::showbase << numberAndFactor[1] << '.' << Timer::endl;

    return 0;
}
