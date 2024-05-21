#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <Aeu.h>

using primeType = unsigned;
std::vector<primeType> loadPrimes(const std::filesystem::path& fromLocation, std::size_t primesCount) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");

//    const auto primesCount = 1024; std::filesystem::file_size(fromLocation) / sizeof(primeType);
    std::ifstream input(fromLocation, std::ios::binary);

    std::vector<primeType> primes (primesCount);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(primeType));

    return primes;
}

using Uns = Aeu<8416>;

void kernel(Uns* const numberAndFactor,
            const primeType* const primes,
            const std::size_t primesCount,
            const std::size_t primesPerIteration) {
    const std::size_t threadId = 65535,//rand() % 65536,
            threadsCount = 65536,
            verboseThreadId = threadId;

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


        Uns product = base;
        for(std::size_t j = startingPrimeIndex; j < endingPrimeIndex; ++j)
            product *= primes[j];

        const auto candidate = Uns::gcd(Uns::powm(a, product, n) - 1, n);
        if(candidate > 1) {
            char buffer [256] {};
            product.getString<10>(buffer, 256, true, false);
            printf("Thread %u. Found factor for power %s, a = %llu.\n", threadId, buffer, a);
            return;
        }

        if(threadId == verboseThreadId)
            printf("Iteration %u is about to complete in all threads.\n", i);
        if(threadId == verboseThreadId)
            printf("Iteration: %u, starting prime: %u, ending prime: %u.\n", i, startingPrimeIndex, endingPrimeIndex);
    }
}

int main() {
    std::vector<std::size_t> primesIndexes = {
            4495744, 26272001, 3661285, 33914092,
            2856137, 96593623, 750612
    };

    std::size_t primesTotal = 100'000'000;
    std::size_t threadsTotal = 256 * 256;
    std::size_t primesPerThreadTotal = primesTotal / threadsTotal;
    std::size_t primesPerIteration = 512;

    for(const auto index: primesIndexes) {
        for(std::size_t threadId = 0; threadId < 256 * 256; ++threadId) {
            const auto threadStart = threadId * primesPerThreadTotal;

            for(std::size_t threadIteration = 0; threadIteration * primesPerIteration < primesPerThreadTotal; ++threadIteration) {
                const std::size_t iterationStart = threadStart + threadIteration * primesPerIteration;

                if(index > iterationStart && index < iterationStart + primesPerIteration) {
                    std::cout << "Index " << index << " is " << index - iterationStart << "'s in thread " << threadId << ", iteration " << threadIteration << std::endl;
                    break;
                }
            }
        }
    }

    kernel(100'000'000, 512);
    //     203'227'136
}