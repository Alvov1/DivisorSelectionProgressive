#include <iostream>
#include <fstream>
#include <filesystem>
#include "Aesi.h"
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





const struct { unsigned x {}; unsigned y {}; } gridDim = { 4, 1 }, blockDim = { 4, 1 }, blockIdx = { 0, 0 }, threadIdx = { 0, 0 };

const unsigned threadId = 3,//blockDim.x * blockIdx.x + threadIdx.x,
                threadsCount = gridDim.x * blockDim.x,
                iterationBorder = 64,
                bStart = /*2 + blockIdx.x, */rand() % 1024,
                bShift = gridDim.x,
                bMax = 2000000000U;
using Uns = Aesi<2048>;


Uns countEBasic(unsigned B, const std::vector<uint64_t>& primes) {
    auto primeUl = primes[0];

    Uns e = 1;

    for (unsigned pi = 0; primeUl < B && pi < primes.size(); ++pi) {
        const auto power = static_cast<unsigned>(log(static_cast<double>(B)) / log(static_cast<double>(primeUl)));
        e *= static_cast<uint64_t>(pow(static_cast<double>(primeUl), static_cast<double>(power)));
        primeUl = primes[pi + 1];
    }

    return e;
}

Uns countEComplex(const std::vector<uint64_t>& primes, const Uns& number) {

    /* First 32 primes in required powers: 3^5 * 5^5 * 7^4 * 11^3 * 13^3 * 17^3....
     * Currently, 247 bits. */
    constexpr Uns eBase = Uns(1) * 2187 * 3125 * 2401 * 1331 * 2197
            * 4913 * 19*19 * 23*23 * 29*29 * 31*31 * 37*37 * 41*41
            * 43 * 47 * 53 * 59 * 61 * 67 * 71 * 73 * 79 * 83 * 89
            * 97 * 101 * 103 * 107 * 109 * 113 * 127 * 131;

    static unsigned lastGroup = 0;

    Uns e = Uns::power2(static_cast<unsigned>(number.bitCount() / 10)) * eBase;

    for(; e.bitCount() < number.bitCount() * 5; ++lastGroup) {
        const unsigned iFirst = 32 + lastGroup * threadsCount + threadId,
            iSecond = lastGroup * threadsCount + threadId + threadsCount / 2 + 32;
        std::cout << primes[iFirst] << " - " << primes[iSecond] << std::endl;
        e *= primes[iFirst]; e *= primes[iSecond];
    }

    return e;
}

void kernel(const std::vector<uint64_t>& primes, std::pair<Uns, Uns>& numberAndFactor) {
    auto& [n, factor] = numberAndFactor;
    const auto checkFactor = [&n, &factor] (const Uns& candidate) {
        if(candidate < 2 || candidate >= n)
            return false;
        factor = candidate; return true;
    };

    uint64_t a = threadId * iterationBorder + 2;
    for (unsigned B = bStart; B < bMax; B += bShift) {
        const Uns e = countEComplex(primes, n);
        if (e == 1) continue;

        for (unsigned it = 0; it < iterationBorder; ++it) {
            if(!factor.isZero())
                return;

            if(checkFactor(Uns::gcd(Uns::powm(a, e, n) - 1, n)))
                return;

            a += threadsCount * iterationBorder;
        }
    }
}

int main(int argc, const char* const* const argv) {
    try {
        if (argc < 3)
            return std::printf("Usage: %s <number> <primes location>", argv[0]);

        std::pair<Uns, Uns> numberAndFactor = { { std::string_view(argv[1]) }, { 0 } };
        Timer::init() << "Factorizing number " << std::hex << std::showbase << numberAndFactor.first <<
            std::dec << " (" << numberAndFactor.first.bitCount() << " bits)." << Timer::endl;

        const std::vector<uint64_t> primes = loadPrimes(argv[2]);
        Timer::out << "Loaded prime table of " << primes.size() << " elements." << Timer::endl;

        kernel(primes, numberAndFactor);
        if (numberAndFactor.second != 0)
            Timer::out << "Found factor: " << std::hex << std::showbase << numberAndFactor.second << Timer::endl;
        else Timer::out << "Factor is not fount" << Timer::endl;

        if(numberAndFactor.first % numberAndFactor.second != 0)
            Timer::out << "Factor is incorrect." << Timer::endl;
    } catch (const std::exception& e) {
        return std::printf("Failed: %s.\n", e.what());
    }
    return 0;
}
