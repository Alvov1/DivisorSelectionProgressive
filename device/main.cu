#include <iostream>
#include <fstream>
#include <filesystem>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "Aesi.h"
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


using Uns = Aesi<1024>;

__global__
void kernel(Uns* const numberAndFactor, const uint64_t* const primes, std::size_t primesCount) {
    const unsigned threadId = blockDim.x * blockIdx.x + threadIdx.x,
            threadsCount = gridDim.x * blockDim.x,
            iterationsBorder = 1024,
            bStart = 2 + blockIdx.x,
            bShift = gridDim.x,
            bMax = 2000000000U;

    const Uns n = numberAndFactor[0]; Uns* const factor = numberAndFactor + 1;
    const auto checkFactor = [&n, &factor] (const Uns& candidate) {
        if(candidate < 2 || candidate >= n)
            return false;
        factor->tryAtomicSet(candidate); return true;
    };
    
    const auto countE = [] (unsigned B, const uint64_t* const primes, std::size_t primesCount) {
        auto primeUl = primes[0];

        Uns e = 1;
        for(unsigned pi = 0; primeUl < B && pi < primesCount; ++pi) {
            const auto power = static_cast<unsigned>(log(static_cast<double>(B) / log(static_cast<double>(primeUl))));
            e *= static_cast<uint64_t>(pow(static_cast<double>(primeUl), static_cast<double>(power)));
            primeUl = primes[pi + 1];
        }

        return e;
    };

    static constexpr auto overflowBoarder = 5;
    static constexpr auto overflowRatio = 0.9;

    uint64_t a = threadId * iterationsBorder + 2, overflowCount = 0;
    for (unsigned B = bStart; B < bMax; B += bShift) {
        const Uns l = countE(B, primes, primesCount);
        if (l == 1) continue;

        for (unsigned it = 0; it < iterationsBorder; ++it) {
            if(!factor->isZero())
                return;

            /* GCD(a^l mod n, n) */
            if(checkFactor(Uns::gcd(Uns::powm(a, l, n) - 1, n)))
                return;

            a += threadsCount * iterationsBorder;
        }
    }
}

int main(int argc, const char* const* const argv) {
    if(argc < 3)
        return std::printf("Usage: %s <number> <primes location>", argv[0]);

    const Uns number = std::string_view(argv[1]);
    thrust::device_vector<Uns> numberAndFactor = { number, Uns { 0 } };
    Timer::init() << "Factorizing number " << std::hex << std::showbase << number << std::dec << " (" << number.bitCount() << " bits)." << Timer::endl;

    const thrust::device_vector<uint64_t> primes = loadPrimes(argv[2]);
    Timer::out << "Loaded prime table of " << primes.size() << " elements." << Timer::endl;

    kernel<<<256, 256>>>(
            thrust::raw_pointer_cast(numberAndFactor.data()),
            thrust::raw_pointer_cast(primes.data()),
            primes.size());

    const auto code = cudaDeviceSynchronize();
    if (code != cudaSuccess)
        return std::printf("Kernel launch failed: %s.\n", cudaGetErrorString(code));
    Timer::out << "Kernel completed. Founded factor: " << std::hex << std::showbase << numberAndFactor[1] << '.' << Timer::endl;

    return 0;
}
