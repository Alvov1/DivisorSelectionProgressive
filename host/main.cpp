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
using TestPrecision = Aeu<ValuePrecision::getBitness() + 32 * 4>;

bool checkValue(const TestPrecision& value) {
    return false; }



void processPrime(const std::vector<primeType>& primes, std::size_t threadId, std::size_t threadsCount) {
    /* Using base with some smallest primes: 2^9 * 3^6 * 5^5 * 7^4 * 11^2 * 13^2 *
     * 17 * 19 * 21 * 23 * 27 * 29 * 31 * 37 = 57267964353600000 (56 bit). */
    constexpr TestPrecision startingBase =
            TestPrecision(512) * 729
            * TestPrecision(3125) * 2401 * 121 * 169
            * TestPrecision(17) * 19 * 21 * 23
            * TestPrecision(27) * 29 * 31 * 37;
    constexpr auto stackSizeTestPrecision = 10,
        primesPositionShift = 1;

    for(std::size_t i = threadId; i < primes.size(); i += threadsCount) {
        /* Storing in stack: prime number with its position in array. */
        PrimeStack<primeType, stackSizeTestPrecision + 2> stack {};
        std::size_t primesPosition = threadId;

        do {
            if (stack.getSize() >= stackSizeTestPrecision) {
                const auto test = stack.unite<TestPrecision>() * startingBase * primes[i];
                if(checkValue(test)) return;
            }

            if (stack.getSize() >= stackSizeTestPrecision || primesPosition >= primes.size()) {
                std::size_t checkedLast = primes.size();
                do {
                    primesPosition = stack.pop();
                    checkedLast -= primesPositionShift;
                } while (primesPosition == checkedLast);

                primesPosition += primesPositionShift;
//                if(primesPosition % threadsCount == threadId)
//                    primesPosition += primesPositionShift;

                stack.push(primes[primesPosition], primesPosition);

                primesPosition += primesPositionShift;
//                if(primesPosition % threadsCount == threadId)
//                    primesPosition += primesPositionShift;
            } else {
                stack.push(primes[primesPosition], primesPosition);

                primesPosition += primesPositionShift;
//                if(primesPosition % threadsCount == threadId)
//                    primesPosition += primesPositionShift;
            }

        } while (!stack.isEmpty());
    }
}

int main() {
//    const std::filesystem::path primesLocation = "../../all-primes-32-bit.bin";
//    const auto primes = loadPrimes(primesLocation, 100);
//    processPrime(primes, 4, 7);

//    constexpr TestPrecision startingBase =
//            TestPrecision(512) * 729
//            * TestPrecision(3125) * 2401 * 121 * 169
//            * TestPrecision(17) * 19 * 21 * 23
//            * TestPrecision(27) * 29 * 31 * 37;
//
//    const Aeu128 a = 3;
//    const Aeu128 n = "0x66ee64b74c3d88a79c57064d8365";
//    const Aeu128 power = Aeu128(2) * 3 * 7 * 83 * 257 * 76840783; //startingBase * 83 * 257 * 76840783;
//
//    //2 × 3 × 7 × 83 × 257 × 76840783
//    //2 × 3 × 7 × 83 × 257 × 76840783
//
//    std::cout << Aeu128::gcd(Aeu128::powm(a, power, n) - 1, n) << std::endl;

//    std::cout << std::is_same_v<uint32_t, unsigned long> << std::endl;
//    std::cout << sizeof(unsigned long) << std::endl;

    std::cout << (1ULL << 32) << std::endl;

    //  11396 22080 // разница 10684
    //  32481 195495 // разница 163014

    std::size_t x = 123;

    for(std::size_t k = 0; k < 999999999; ++k) {
        for (unsigned a = 0; a < 999999999; ++a) {
            for (unsigned b = 0; b < 999999999; ++b) {
                if(k + a * x == 10684 && k + b * x == 163014)
                    std::cout << "k = " << k << std::endl << "a = " << a << std::endl << "b = " << b << std::endl;
            }
        }
    }
}