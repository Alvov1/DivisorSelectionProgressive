#include <iostream>
#include <fstream>
#include <filesystem>
#include <bitset>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cuda/std/array>
#include "Aeu.h"
#include "Timer.h"

using Uns = Aeu<4096>;

using primeType = uint64_t;
std::size_t binary_search_find_index(const thrust::host_vector<primeType>& v, primeType data) {
    auto it = std::lower_bound(v.begin(), v.end(), data);
    if (it == v.end()) {
        return 0;
    } else {
        return std::distance(v.begin(), it);
    }
}


thrust::host_vector<primeType> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");


    const auto primesCount = std::filesystem::file_size(fromLocation) / sizeof(primeType);
    std::cout << "Loading prime table of " << primesCount << " elements ("
        << std::filesystem::file_size(fromLocation) << " bytes)" << std::endl;
    std::ifstream input(fromLocation, std::ios::binary);

    thrust::host_vector<primeType> primes (primesCount);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(primeType));

    return primes;
}

constexpr auto getNumbersBitness = [] (primeType number) {
    return sizeof(primeType) * 8 - std::countl_zero(number);
};

__global__
void kernel(Uns* const numberFactorBase,
            const primeType* const primes,
            const unsigned primesCount,
            const unsigned primesPerIteration,
            const unsigned verboseThreadId) {
    const unsigned threadId = blockDim.x * blockIdx.x + threadIdx.x,
            threadsCount = gridDim.x * blockDim.x,
            logicalThreadId = (threadId + verboseThreadId) % threadsCount;

    const Uns& n = numberFactorBase[0], a = 2, &base = numberFactorBase[2];
    Uns& factor = numberFactorBase[1];

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
            printf("Thread %u. Found factor.\n", logicalThreadId);
            factor.tryAtomicSet(candidate);
            return;
        }

        if(logicalThreadId == verboseThreadId)
            printf("Starting prime: %u, ending prime: %u. \t---\t Iteration %u complete\n", startingPrimeIndex, endingPrimeIndex, i);
    }

    if(logicalThreadId == verboseThreadId)
        std::printf("Verbose thread exited.\n");
}



int showUsage(const char* const exe) {
    return std::printf("Usage: %s <number> <factorization base> <primes location> <threads> <verbose thread id>\n"
                       "\tNumber - number to factorize\n"
                       "\tFactorization base - product of the smallest primes. IT MUST CONTAIN ALL FACTORS FOR THE (p-1) EXCEPT FOR THE BIGGEST ONE\n"
                       "\tPrimes location - file with all available primes, SORTED IN ASCENDING ORDER\n"
                       "\tThreads - amount of threads (used both for amount of blocks and amount of threads in one block)\n"
                       "\tVerbose thread id - id of the only one thread which writes logs in the stdout\n\n", exe);
}




int main(int argc, const char* const* const argv) {
    try {
        if(argc < 6)
            return showUsage(argv[0]);

        /* 1. Load number for factorize and factorization base. */
        const Uns number = std::string_view(argv[1]), factorBase = std::string_view(argv[2]);
        thrust::device_vector<Uns> numberFactorBase = thrust::host_vector{ number, Uns { 0u }, factorBase };
        Timer::init() << "Factorizing number 0x" << std::hex << number
                      << ". Using factor base 0x" << factorBase << std::dec << Timer::endl;

        /* 2. Load prime table and filter with primes of only required bitness. */
        const std::filesystem::path primesLocation (argv[3]);
        const thrust::device_vector<primeType> primes = loadPrimes(primesLocation);
        Timer::out << "Loaded prime table, which was filtered for factors " << remainingFactorBitness
                   << " bit length only: " << primes.size() << " elements." << Timer::endl;
        Timer::out << "First prime: " << primes[0] << " (" << std::bitset<sizeof(primeType) * 8>(primes[0])
                   << "), last prime: " << primes.back() << " (" << std::bitset<sizeof(primeType) * 8>(primes.back()) << ")." << Timer::endl;

        /* 3. Calculate the number of primes for selection for each thread. */
        const std::size_t threads = std::stoul(argv[4]), verboseThreadId = std::stoul(argv[5]);
        const unsigned primesPerIteration = (Uns::getBitness() - factorBase.bitCount() - 32) / 32;
        if(primesPerIteration < 50)
            Timer::out << "Warning: selection window may be too small for enough efficiency (taking " << primesPerIteration
                       << " primes at each iteration). Recompilation with greater precision length may be required." << Timer::endl;

        /* 4. Launch kernel with all the data. */
        const auto timePoint = std::chrono::system_clock::now();
        const auto timeT = std::chrono::system_clock::to_time_t(timePoint);
        Timer::out << std::ctime(&timeT) << " Starting kernel <<<" << threads << ", " << threads << ">>>. Using bitness "
                   << Uns::getBitness() << ". Primes per iteration: \n" << primesPerIteration << Timer::endl;

        /* 5. Wait until kernel completes. */
        kernel<<<threads, threads>>>(
                thrust::raw_pointer_cast(numberFactorBase.data()),
                thrust::raw_pointer_cast(primes.data()),
                primes.size(),
                primesPerIteration,
                verboseThreadId);
        const auto code = cudaDeviceSynchronize();
        if(code != cudaSuccess)
            return std::printf("Kernel launch failed: %s.\n", cudaGetErrorString(code));


        /* 6. Review the results. */
        const Uns factor = numberFactorBase[1];
        if(factor != 0u && number % factor == 0u)
            Timer::out << "Kernel completed. Founded factor: 0x" << std::hex << factor << '.' << Timer::endl;
        else Timer::out << "Kernel completed, but factor 0x" << std::hex << factor << " is incorrect." << Timer::endl;

    } catch (const std::exception& e) {
        return std::printf("Failed: %s.\n", e.what());
    }
}