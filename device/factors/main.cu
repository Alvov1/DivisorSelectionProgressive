#include <iostream>
#include <fstream>
#include <filesystem>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cuda/std/array>
#include "Aeu.h"
#include "Timer.h"

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








/* 128 primes + 192 bit base + 32 bit additional -> 4320 */
using Uns = Aeu<2272>;

__global__
void kernel(Uns* const numberAndFactor,
            const primeType* const primes,
            const unsigned primesCount,
            const unsigned primesPerIteration,
            const unsigned verboseThreadId) {
    const unsigned threadId = blockDim.x * blockIdx.x + threadIdx.x,
            threadsCount = gridDim.x * blockDim.x,
            logicalThreadId = (threadId + verboseThreadId) % threadsCount;

    const Uns& n = numberAndFactor[0], a = 2;
    Uns& factor = numberAndFactor[1];

    Uns base = Uns(1) * 512 * 729 * 3125 * 2401 * 121 *
               169 * 17 * 19 * 21 * 23 * 27 * 29
               * 31 * 37 * 83 * 131 * 199 * 257
               * 317 * 1051 * 1297 * 1637 * 1933;

    /* Количество простых чисел в области потока. */
    const unsigned primesPerThread = primesCount / threadsCount;

    /* Количество итераций, которое необходимо чтобы перебрать все эти числа. */
    const unsigned iterationsAmount = primesPerThread / primesPerIteration
                                      + (primesPerThread % primesPerIteration == 0 ? 0 : 1);

    /* Индекс первого числа (начала области потока). */
    const unsigned firstPrimeInThread = logicalThreadId * primesPerThread;

    if(logicalThreadId == verboseThreadId) {
        printf("\n\nVerbose thread: %u\n", logicalThreadId);
        printf("Factorizing number with bitness %u\n", n.bitCount());
        printf("Base bitness %u\n", base.bitCount());
        printf("Primes per thread: %u\n", primesPerThread);
        printf("Primes per iteration: %u\n", primesPerIteration);
        printf("Iterations amount: %u\n", iterationsAmount);
        printf("First prime in thread: %u\n\n\n", firstPrimeInThread);
    }

    for(unsigned i = 0; i < iterationsAmount; ++i) {
        /* Индекс первого числа в этой итерации. */
        const unsigned startingPrimeIndex = firstPrimeInThread + (i * primesPerIteration);

        /* Индекс последнего числа в этой итерации. */
        const unsigned endingPrimeIndex = (startingPrimeIndex + primesPerIteration < firstPrimeInThread + primesPerThread ?
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
            printf("Thread %u. Found factor.\n", logicalThreadId);

            factor.tryAtomicSet(candidate);
            return;
        }

        if(logicalThreadId == verboseThreadId)
            printf("Starting prime: %u, ending prime: %u. \t---\t Iteration %u complete\n", startingPrimeIndex, endingPrimeIndex, i);
    }
}











int main(int argc, const char* const* const argv) {
    try {
        if(argc < 5)
            return std::printf("Usage: %s <number> <primes location> <threads> <verbose thread id>", argv[0]);

        const Uns number = std::string_view(argv[1]);
        thrust::device_vector<Uns> numberAndFactor = thrust::host_vector{ number, Uns { 0u } };
        Timer::init() << "Factorizing number 0x" << std::hex << number << '.' << std::dec << Timer::endl;




        const thrust::device_vector<primeType> primes = loadPrimes(argv[2]);
        Timer::out << "Loaded table of primes with " << primes.size() << " elements." << Timer::endl;





        const std::size_t threads = std::stoul(argv[3]);
        const unsigned primesPerIteration = (Uns::getBitness() - 224) / 32;
        const unsigned verboseThreadId = std::stoul(argv[4]);

        const auto timePoint = std::chrono::system_clock::now();
        const auto timeT = std::chrono::system_clock::to_time_t(timePoint);
        Timer::out << std::ctime(&timeT) << " Starting kernel <<<" << threads << ", " << threads << ">>>. Using bitness "
                   << Uns::getBitness() << ". Primes per iteration: \n" << primesPerIteration << Timer::endl;





        kernel<<<threads, threads>>>(
                thrust::raw_pointer_cast(numberAndFactor.data()),
                thrust::raw_pointer_cast(primes.data()),
                primes.size(),
                primesPerIteration,
                verboseThreadId);






        const auto code = cudaDeviceSynchronize();
        if(code != cudaSuccess)
            return std::printf("Kernel launch failed: %s.\n", cudaGetErrorString(code));





        const Uns factor = numberAndFactor[1];
        if(factor != 0u && number % factor == 0u)
            Timer::out << "Kernel completed. Founded factor: 0x" << std::hex << factor << '.' << Timer::endl;
        else Timer::out << "Kernel completed, but factor 0x" << std::hex << factor << " is incorrect." << Timer::endl;

    } catch (const std::exception& e) {
        return std::printf("Failed: %s.\n", e.what());
    }
}