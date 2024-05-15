#ifndef DIVISORSELECTIONPROGRESSIVE_SElECTION_KERNEL_H
#define DIVISORSELECTIONPROGRESSIVE_SElECTION_KERNEL_H

#include "../../host/simple-circular-buffer.h"

using UCalculations = Aeu<2048>;
using UBase = Aeu<256>;
using primeType = unsigned;

__global__
void kernel(UBase* const numberAndFactor, const unsigned* const primes, std::size_t primesCount, std::size_t iterationsPerPrime) {
    const unsigned threadId = blockDim.x * blockIdx.x + threadIdx.x,
            threadsCount = gridDim.x * blockDim.x,
            bStart = 2 + blockIdx.x,
            bShift = gridDim.x,
            bMax = 2000000000U;

    const Uns n = numberAndFactor[0]; Uns* const factor = numberAndFactor + 1;
    const auto checkFactor = [&n, &factor] (const UCalculations& candidate) {
        if(candidate < 2 || candidate >= n)
            return false;
        factor->tryAtomicSet(candidate); return true;
    };
    constexpr auto getNumbersBitness = [] (primeType prime) {
        return (sizeof(primeType) * 8 - __clz(prime));
    };

    /* Average amount of prime number selected for the thread. */
    std::size_t numbersPerThread = primesCount / threadsCount;

    /* Power base, containing the smallest primes in some powers. */
    constexpr UCalculations base = UCalculations::power2(9) * 2187 * 3125 * 2401; // 2^9 * 3^7 * 5^5 * 7^4
    const UCalculations a = 3;

    /* Approximate number of primes required to overflow testing bitness
     * Ex: testing bitness 2048 -> 86 numbers of bitness 24. */
    constexpr auto approximateNumberOfPrimesPerTestingBitness = static_cast<std::size_t>(UCalculations::getBitness() / (sizeof(primeType) * 4));

    /* Shift in prime table according to which we're taking primes. */
    const auto primesShift = static_cast<std::size_t>(primesCount / 60);

    for(std::size_t i = 0; i < numbersPerThread; ++i) {
        const std::size_t currentNumberIndex = (threadId + i * numbersPerThread) % primesCount;

        CircularBuffer<primeType, approximateNumberOfPrimesPerTestingBitness> buffer {};
        UCalculations e = base;

        std::size_t nextPrimeIndex = threadId % primesShift;

        for(std::size_t j = 0; j < iterationsPerPrime; ++j) {
            const primeType nextPrime = primes[nextPrimeIndex];
            const auto nextPrimeBitness = getNumbersBitness(nextPrime);

            /* If in the number there is enough space for new prime. */
            const auto currentBitness = e.bitCount();
            if(currentBitness + nextPrimeBitness + 10 < UCalculations::getBitness()) {
                e *= nextPrime;
                buffer.push_back(nextPrime);
                nextPrimeIndex = (nextPrimeIndex + primesShift) % primesCount;
            } else {
                /* If the number's overflowed. */
                /* 1. Check. */
                if(checkFactor(UBase::gcd(UBase::powm(a, e, n) - 1, n))) {
                    char buffer [512] {};
                    e.getString<10>(buffer, 512, true, false);
                    printf("Thread %u. Found factor for power %s, a = %llu.\n", threadId, buffer, a);
                    return;
                }


                /* 2. Erase bitness. */
                for(std::size_t t = 0; e.bitCount() + nextPrimeBitness + 10 >= UCalculations::getBitness() && t < approximateNumberOfPrimesPerTestingBitness; ++t) {
                    e /= UCalculations(buffer.pop_front());
                }
            }
        }
    }
}

#endif //DIVISORSELECTIONPROGRESSIVE_SElECTION_KERNEL_H
