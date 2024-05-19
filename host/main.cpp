#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <Aeu.h>

#include "PrimeStack.h"

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

//void kernel(const primeType* const primes, std::size_t primesCount,
//            std::size_t threadId, std::size_t threadsCount) {
//    for(std::size_t i = threadId; i < primesCount; i += threadsCount) {
//        const auto prime = primes[i];
//        processPrime(primes, primesCount, i, 30);
//    }
//}

using ValuePrecision = Aeu128;
using TestPrecision = Aeu<ValuePrecision::getBitness() * 2 + 32>;

bool checkValue(const TestPrecision& value) {
    return false; }



void processPrime(const std::vector<primeType>& primes, std::size_t threadId, std::size_t threadsCount) {
    /* Using base with some smallest primes: 2^9 * 3^6 * 5^5 * * 11^2 * 13^2 *
     * 17 * 19 * 21 = 57267964353600000 (56 bit). */
    constexpr TestPrecision startingBase = TestPrecision(512) * 729 * TestPrecision(3125) * 2401 * 121 * 169;
    constexpr auto stackSizeTestPrecision = TestPrecision::getBitness() / 8;

    for(std::size_t i = threadId; i < primes.size(); i += threadsCount) {
        /* Storing in stack: prime number with its position in array. */
        PrimeStack<primeType, stackSizeTestPrecision + 2> stack {};
        std::size_t primesPosition = 0;

        do {
            if (stack.getSize() >= stackSizeTestPrecision) {
                const auto test = stack.unite<TestPrecision>() * startingBase * primes[i];
                if(checkValue(test)) return;
            }

            if (stack.getSize() >= stackSizeTestPrecision || primesPosition == primes.size()) {
                std::size_t checkedLast = primes.size() - 1;
                do {
                    primesPosition = stack.pop();
                } while (primesPosition == checkedLast--);

                if(++primesPosition % threadsCount == threadId)
                    ++primesPosition;
                stack.push(primes[primesPosition], primesPosition);
                if(++primesPosition % threadsCount == threadId)
                    ++primesPosition;
            } else {
                stack.push(primes[primesPosition], primesPosition);
                if(++primesPosition % threadsCount == threadId)
                    ++primesPosition;
            }

        } while (!stack.isEmpty());
    }
}

int main() {
    const std::filesystem::path primesLocation = "../../all-primes-32-bit.bin";
    const auto primes = loadPrimes(primesLocation, 1024);

    processPrime(primes, 4, 7);
}