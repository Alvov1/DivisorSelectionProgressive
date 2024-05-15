#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include "Aeu.h"
#include "Timer.h"
#include "simple-circular-buffer.h"

const struct { unsigned x {}; unsigned y {}; } gridDim = { 256, 1 }, blockDim = { 256, 1 }, blockIdx = { 0, 0 }, threadIdx = { 0, 0 };

const unsigned threadId = 0,//3,//blockDim.x * blockIdx.x + threadIdx.x,
threadsCount = gridDim.x * blockDim.x,
        iterationBorder = 64,
        bStart = 2 + blockIdx.x,
        bShift = gridDim.x,
        bMax = 2'000'000'000U;
using Uns = Aeu512;

using primeType = unsigned;
std::vector<primeType> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");

    const auto primesCount = 3388003;//std::filesystem::file_size(fromLocation) / sizeof(primeType);
    std::ifstream input(fromLocation, std::ios::binary);

    std::vector<primeType> primes (primesCount);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(primeType));

    return primes;
}

void kernel(const std::vector<primeType>& primes, std::pair<Uns, Uns>& numberAndFactor) {
    auto& [n, factor] = numberAndFactor;
    const auto checkFactor = [&n, &factor] (const Uns& candidate) {
        if(candidate < 2u || candidate >= n)
            return false;
        factor = candidate; return true;
    };
    constexpr auto getNumbersBitness = [] (primeType prime) {
#ifndef __CUDACC__
        return (sizeof(primeType) * 8 - std::countl_zero(prime));
#else
        return (sizeof(primeType) * 8 - __clz(lastBlock));
#endif
    };
    const auto primesCount = primes.size();
    const std::size_t iterationsPerPrime = 2048;

    /* Average amount of prime number selected for the thread. */
    std::size_t numbersPerThread = primesCount / threadsCount;

    /* Power base, containing the smallest primes in some powers. */
    constexpr Uns base = Uns::power2(9) * 2187 * 3125 * 2401, a = 3; // 2^9 * 3^7 * 5^5 * 7^4
    constexpr auto baseBitness = base.bitCount();

    /* Approximate number of primes required to overflow testing bitness
     * Ex: testing bitness 2048 -> 86 numbers of bitness 24. */
    constexpr auto approximateNumberOfPrimesPerTestingBitness = static_cast<std::size_t>(Uns::getBitness() / (sizeof(primeType) * 4));

    /* Shift in prime table according to which we're taking primes. */
    const auto primesShift = static_cast<std::size_t>(primesCount / 30);

    for(std::size_t i = 0; i < numbersPerThread; ++i) {
        const std::size_t currentNumberIndex = (threadId + i * numbersPerThread) % primesCount;

        CircularBuffer<primeType, approximateNumberOfPrimesPerTestingBitness> buffer {};
        Uns e = base;

        std::size_t nextPrimeIndex = threadId % primesShift;

        for(std::size_t j = 0; j < iterationsPerPrime; ++j) {
            const primeType nextPrime = primes[nextPrimeIndex];
            const auto nextPrimeBitness = getNumbersBitness(nextPrime);

            /* If in the number there is enough space for new prime. */
            const auto currentBitness = e.bitCount();
            if(currentBitness + nextPrimeBitness + 10 < Uns::getBitness()) {
                e *= nextPrime;
                buffer.push_back(nextPrime);
                nextPrimeIndex = (nextPrimeIndex + primesShift) % primesCount;
            } else {
                /* If the number's overflowed. */
                /* 1. Check. */
                buffer.introspect();

                if(checkFactor(Uns::gcd(Uns::powm(a, e, n) - 1, n))) {
                    std::cout << "Thread " << threadId << ". Found factor for power " << e << ", a = " << a << std::endl;
                    return;
                }


                /* 2. Erase bitness. */
                for(std::size_t t = 0; e.bitCount() + nextPrimeBitness + 10 >= Uns::getBitness() && t < approximateNumberOfPrimesPerTestingBitness; ++t) {
                    e /= Uns(buffer.pop_front());
                }
            }
        }
    }
}

int main() {
    const std::filesystem::path primes = "../../all-primes-32-bit.bin";
    const Uns tNumber = "0xb8dfad5e20f1c2d";

    std::pair<Uns, Uns> numberAndFactor = { tNumber, { 0 } };
    Timer::init() << "Factorizing number " << std::hex << std::showbase << numberAndFactor.first <<
                  std::dec << " (" << numberAndFactor.first.bitCount() << " bits)." << Timer::endl;

    const auto values = loadPrimes(primes);
    Timer::out << "Loaded table of smooth with " << values.size() << " elements." << Timer::endl;

    kernel(values, numberAndFactor);
    if (numberAndFactor.second != 0 && numberAndFactor.first % numberAndFactor.second == 0)
        Timer::out << "Found factor: " << std::hex << std::showbase << numberAndFactor.second << Timer::endl;
    else Timer::out << "Factor is not fount or incorrect" << Timer::endl;
    return 0;
}