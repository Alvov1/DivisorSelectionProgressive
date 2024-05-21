#include <iostream>
#include <fstream>
#include <filesystem>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "Aeu.h"
#include "Timer.h"
#include "PrimeStack.h"

using primeType = unsigned;
thrust::host_vector<primeType> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");

    const auto primesCount = (std::filesystem::file_size(fromLocation) / sizeof(primeType)) / 2;
    std::ifstream input(fromLocation, std::ios::binary);

    thrust::host_vector<primeType> primes (primesCount);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(primeType));

    return primes;
}








/* 256 primes + 192 bit base + 32 bit additional */
using Uns = Aeu<8416>;

__global__
void kernel(Uns* const numberAndFactor,
            const primeType* const primes,
            const std::size_t primesCount,
            const std::size_t primesPerIteration) {
    const unsigned threadId = blockDim.x * blockIdx.x + threadIdx.x,
            threadsCount = gridDim.x * blockDim.x,
            verboseThreadId = 0;

    const Uns& n = numberAndFactor[0], a = 2;
    Uns& factor = numberAndFactor[1];
    const Uns base = Uns("0x104e822c30254e465a6cb92d99617bbca06b53f5d200"); //(173 bits)

    /* Количество простых чисел в области потока. */
    const std::size_t primesPerThread = primesCount / threadsCount;
    if(threadId == verboseThreadId)
        printf("Primes per thread: %u\n", primesPerThread);

    /* Количество итераций, которое необходимо чтобы перебрать все эти числа. */
    const std::size_t iterationsAmount = primesPerThread / primesPerIteration
                                         + (primesPerThread % primesPerIteration == 0 ? 0 : 1);
    if(threadId == verboseThreadId)
        printf("Iterations amount: %u\n", iterationsAmount);

    /* Индекс первого числа (начала области потока). */
    const std::size_t firstPrimeInThread = threadId * primesPerThread;
    if(threadId == verboseThreadId)
        printf("First prime in thread: %u\n", firstPrimeInThread);

    for(std::size_t i = 0; i < iterationsAmount; ++i) {
        /* Индекс первого числа в этой итерации. */
        const std::size_t startingPrimeIndex = firstPrimeInThread + i * primesPerIteration;

        /* Индекс последнего числа в этой итерации. */
        const std::size_t endingPrimeIndex = (startingPrimeIndex + primesPerIteration < firstPrimeInThread + primesPerThread ?
                                              startingPrimeIndex + primesPerIteration
                                                                                                                             :
                                              firstPrimeInThread + primesPerThread);

        if(!factor.isZero())
            return;

        Uns product = base;
        for(std::size_t j = startingPrimeIndex; j < endingPrimeIndex; ++j)
            product *= primes[j];

        const auto candidate = Uns::gcd(Uns::powm(a, product, n) - 1, n);
        if(candidate > 1) {
            char buffer [256] {};
            product.getString<10>(buffer, 256, true, false);
            printf("Thread %u. Found factor for power %s, a = %llu.\n", threadId, buffer, a);

            factor.tryAtomicSet(candidate);
            return;
        }

        if(threadId == verboseThreadId)
            printf("Iteration %u is about to complete in all threads.\n", i);
        if(threadId == verboseThreadId)
            printf("Iteration: %u, starting prime: %u, ending prime: %u.\n", i, startingPrimeIndex, endingPrimeIndex);
    }
}











int main(int argc, const char* const* const argv) {
    try {
        if(argc < 5)
            return std::printf("Usage: %s <number> <primes location> <threads> <threads per iteration>", argv[0]);

        const ValuePrecision number = std::string_view(argv[1]);
        thrust::device_vector<ValuePrecision> numberAndFactor = thrust::host_vector{ number, ValuePrecision { 0u } };
        Timer::init() << "Factorizing number 0x" << std::hex << number << '.' << std::dec << Timer::endl;




        const thrust::device_vector<primeType> primes = loadPrimes(argv[2]);
        Timer::out << "Loaded table of primes with " << primes.size() << " elements." << Timer::endl;





        const auto threads = std::stoul(argv[3]), iterations = std::stoul(argv[4]);
        const auto timePoint = std::chrono::system_clock::now();
        const auto timeT = std::chrono::system_clock::to_time_t(timePoint);
        Timer::out << std::ctime(&timeT) << " Starting kernel <<<" << threads << ", " << threads << ">>>. Using bitness "
                   << ValuePrecision::getBitness() << ". Primes per iteration: " << iterations << Timer::endl;





        kernel<<<threads, threads>>>(
                thrust::raw_pointer_cast(numberAndFactor.data()),
                thrust::raw_pointer_cast(primes.data()),
                primes.size(),
                iterations);






        const auto code = cudaDeviceSynchronize();
        if(code != cudaSuccess)
            return std::printf("Kernel launch failed: %s.\n", cudaGetErrorString(code));





        const ValuePrecision factor = numberAndFactor[1];
        if(factor != 0u && number % factor == 0u)
            Timer::out << "Kernel completed. Founded factor: 0x" << std::hex << number << '.' << Timer::endl;
        else Timer::out << "Kernel completed, but factor 0x" << std::hex << factor << " is incorrect." << Timer::endl;

    } catch (const std::exception& e) {
        return std::printf("Failed: %s.\n", e.what());
    }
}