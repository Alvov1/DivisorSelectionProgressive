#include <iostream>
#include <fstream>
#include <filesystem>
#include "AesiMultiprecision.h"
#include "Timer.h"

std::vector<uint64_t> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");

    std::ifstream input(fromLocation, std::ios::binary);
    uint64_t buffer {}; input.read(reinterpret_cast<char*>(&buffer), sizeof(uint64_t));

    std::vector<uint64_t> primes (buffer);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(uint64_t));

    return primes;
}

constexpr auto blocksNumber = 64, threadsPerBlock = 64;
const struct { unsigned x {}; unsigned y {}; } gridDim = { 64, 1 }, blockDim = { 64, 1 }, blockIdx = { 0, 0 }, threadIdx = { 0, 0 };

const unsigned threadId = blockDim.x * blockIdx.x + threadIdx.x,
                threadsCount = gridDim.x * blockDim.x,
                iterationsBorder = 1024,
                bStart = 2 + blockIdx.x,
                bShift = gridDim.x,
                bMax = 2000000000U;

void kernel(const std::vector<uint64_t>& primes, std::pair<Aesi<512>, Aesi<512>>& numberAndFactor) {
    auto& [n, factor] = numberAndFactor;
    const auto checkFactor = [&n, &factor] (const Aesi<512>& candidate) {
        if(candidate < 2 || candidate >= n)
            return false;
        factor = candidate; return true;
    };

    Aesi a = threadId * iterationsBorder + 2;
    for (unsigned B = bStart; B < bMax; B += bShift) {
        auto primeUl = primes[0];

        Aesi e = 1;
        for (unsigned pi = 0; primeUl < B; ++pi) {
            if(!factor.isZero()) return;
            const auto power = static_cast<unsigned>(log(static_cast<double>(B)) / log(static_cast<double>(primeUl)));
            e *= static_cast<uint64_t>(pow(static_cast<double>(primeUl), static_cast<double>(power)));
            primeUl = primes[pi + 1];
        }

        if (e == 1) continue;

        for (unsigned it = 0; it < iterationsBorder; ++it) {
            if(!factor.isZero())
                return;

            if(checkFactor(Aesi<512>::gcd(a, n)))
                return;

            if(checkFactor(Aesi<512>::gcd(Aesi<512>::powm(a, e, n) - 1, n)))
                return;

            a += threadsCount * iterationsBorder;
        }
    }
}

int main(int argc, const char* const* const argv) {
    try {
        if (argc < 3)
            return std::printf("Usage: %s <number> <primes location>", argv[0]);

        std::pair<Aesi<512>, Aesi<512>> numberAndFactor = { { std::string_view(argv[1]) }, { 0 } };
        Timer::init() << "Factorizing number " << std::hex << std::showbase << numberAndFactor.first << std::dec << '.' << Timer::endl;

        const std::vector<uint64_t> primes = loadPrimes(argv[2]);
        Timer::out << "Loaded prime table of " << primes.size() << " elements." << Timer::endl;

        kernel(primes, numberAndFactor);
        if (numberAndFactor.second != 0)
            Timer::out << "*Kernel* completed. Fount factor: " << std::hex << std::showbase << numberAndFactor.second << Timer::endl;
        Timer::out << "*Kernel* completed. Factor is not fount" << Timer::endl;
    } catch (const std::exception& e) {
        return std::printf("Failed: %s.\n", e.what());
    }
    return 0;
}
