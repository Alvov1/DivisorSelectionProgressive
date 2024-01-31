#include <iostream>
#include <fstream>
#include <filesystem>
#include "Aesi.h"
#include "Timer.h"

const struct { unsigned x {}; unsigned y {}; } gridDim = { 64, 1 }, blockDim = { 64, 1 }, blockIdx = { 0, 0 }, threadIdx = { 0, 0 };

const unsigned threadId = 0,//3,//blockDim.x * blockIdx.x + threadIdx.x,
threadsCount = gridDim.x * blockDim.x,
        iterationBorder = 1024,
        bStart = 2 + blockIdx.x,
        bShift = gridDim.x,
        bMax = 2'000'000'000U;
using Uns = Aesi<512>;

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

Uns countE(unsigned B, const std::vector<uint64_t>& primes) {
    auto primeUl = primes[0];

    Uns e = 1;

    for (unsigned pi = 0; primeUl < B && pi < primes.size(); ++pi) {
        const auto power = static_cast<unsigned>(log(static_cast<double>(B)) / log(static_cast<double>(primeUl)));
        e *= static_cast<uint64_t>(pow(static_cast<double>(primeUl), static_cast<double>(power)));
        primeUl = primes[pi + 1];
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
        const Uns e = countE(B, primes);
        if (e == 1) continue;

        for (unsigned it = 0; it < iterationBorder; ++it) {
            if(!factor.isZero())
                return;

            if(checkFactor(Uns::gcd(Uns::powm(a, e, n) - 1, n)))
                return;

            a += threadsCount * iterationBorder;
        }
    }

    std::cout << "Here :D" << std::endl;
}

int main() {
    const std::filesystem::path primes = "../../primes.txt";
    const Uns tNumber = "0x1ba4b22fb8d11f51de10f84d";

    std::pair<Uns, Uns> numberAndFactor = { tNumber, { 0 } };
    Timer::init() << "Factorizing number " << std::hex << std::showbase << numberAndFactor.first <<
                  std::dec << " (" << numberAndFactor.first.bitCount() << " bits)." << Timer::endl;

    const auto values = loadPrimes(primes);
    Timer::out << "Loaded table of primes with " << values.size() << " elements." << Timer::endl;

    kernel(values, numberAndFactor);
    if (numberAndFactor.second != 0 && numberAndFactor.first % numberAndFactor.second == 0)
        Timer::out << "Found factor: " << std::hex << std::showbase << numberAndFactor.second << Timer::endl;
    else Timer::out << "Factor is not fount or incorrect" << Timer::endl;
    return 0;
}