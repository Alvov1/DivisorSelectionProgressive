#ifndef DIVISORSELECTIONPROGRESSIVE_POWERS_KERNEL_H
#define DIVISORSELECTIONPROGRESSIVE_POWERS_KERNEL_H

template<std::size_t bitness>
__global__
void kernel(Aesi<bitness>* const numberAndFactor, const Aesi<bitness>* const powers, std::size_t powersCount, std::size_t iterationsBorder) {
    const unsigned threadId = blockDim.x * blockIdx.x + threadIdx.x,
            threadsCount = gridDim.x * blockDim.x;
    using Uns = Aesi<bitness>;

    const Uns n = numberAndFactor[0]; Uns* const factor = numberAndFactor + 1;
    const auto checkFactor = [&n, &factor] (const Uns& candidate) {
        if(candidate < 2 || candidate >= n)
            return false;
        factor->tryAtomicSet(candidate); return true;
    };

    uint64_t a = threadId * iterationsBorder + 2;
    for (unsigned i = threadId; i < powersCount; i += blockDim.x) {
        for (unsigned it = 0; it < iterationsBorder; ++it) {
            if(!factor->isZero())
                return;

            if(checkFactor(Uns::gcd(Uns::powm(a, powers[i], n) - 1, n)))
                return;

            a += threadsCount * iterationsBorder;
        }
    }
}

#endif //DIVISORSELECTIONPROGRESSIVE_POWERS_KERNEL_H
