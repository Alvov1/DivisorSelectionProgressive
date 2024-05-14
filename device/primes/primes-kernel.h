#ifndef DIVISORSELECTIONPROGRESSIVE_PRIMES_KERNEL_H
#define DIVISORSELECTIONPROGRESSIVE_PRIMES_KERNEL_H

using Uns = Aeu<512>;

__device__
        Uns countE(unsigned B, const unsigned* const primes, std::size_t primesCount) {
    auto primeUl = primes[0];

    Uns e = 1;

    for (unsigned pi = 0; primeUl < B && pi < primesCount; ++pi) {
        const auto power = static_cast<unsigned>(log(static_cast<double>(B)) / log(static_cast<double>(primeUl)));
        e *= static_cast<unsigned>(pow(static_cast<double>(primeUl), static_cast<double>(power)));
        primeUl = primes[pi + 1];
    }

    return e;
}

__global__
void kernel(Uns* const numberAndFactor, const unsigned* const primes, std::size_t primesCount, std::size_t iterationsBorder) {
    const unsigned threadId = blockDim.x * blockIdx.x + threadIdx.x,
            threadsCount = gridDim.x * blockDim.x,
            bStart = 2 + blockIdx.x,
            bShift = gridDim.x,
            bMax = 2000000000U;

    const Uns n = numberAndFactor[0]; Uns* const factor = numberAndFactor + 1;
    const auto checkFactor = [&n, &factor] (const Uns& candidate) {
        if(candidate < 2 || candidate >= n)
            return false;
        factor->tryAtomicSet(candidate); return true;
    };

    uint64_t a = threadId * iterationsBorder + 2;
    for (unsigned B = bStart; B < bMax; B += bShift) {
        const Uns e = countE(B, primes, primesCount);
        if (e == 1) continue;

        for (unsigned it = 0; it < iterationsBorder; ++it) {
            if(!factor->isZero())
                return;

            if(checkFactor(Uns::gcd(Uns::powm(a, e, n) - 1, n))) {
                char buffer [256] {};
                e.getString<10>(buffer, 256, true, false);
                printf("Thread %u. Found factor for power %s, a = %llu.\n", threadId, buffer, a);
                return;
            }

            a += threadsCount * iterationsBorder;
        }
    }
}

#endif //DIVISORSELECTIONPROGRESSIVE_PRIMES_KERNEL_H
