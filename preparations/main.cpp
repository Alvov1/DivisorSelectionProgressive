#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <bitset>
#include <Aeu.h>


using Uns = Aeu<8192>;

void kernel(//Uns* const numberAndFactor,
            //const primeType* const primes,
            const std::size_t primesCount,
            const std::size_t primesPerIteration) {
    const std::size_t threadId = 10,//rand() % 65536,
            threadsCount = 65536,
            verboseThreadId = threadId;

//    const Uns& n = numberAndFactor[0], a = 2;
//    Uns& factor = numberAndFactor[1];
//    const Uns base = Uns("0x104e822c30254e465a6cb92d99617bbca06b53f5d200"); //(173 bits)

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


//        Uns product = base;
//        for(std::size_t j = startingPrimeIndex; j < endingPrimeIndex; ++j)
//            product *= primes[j];
//
//        const auto candidate = Uns::gcd(Uns::powm(a, product, n) - 1, n);
//        if(candidate > 1) {
//            char buffer [256] {};
//            product.getString<10>(buffer, 256, true, false);
//            printf("Thread %u. Found factor for power %s, a = %llu.\n", threadId, buffer, a);
//            return;
//        }

//        if(threadId == verboseThreadId)
//            printf("Iteration %u is about to complete in all threads.\n", i);
        if(threadId == verboseThreadId)
            printf("Iteration: %u, starting prime: %u, ending prime: %u.\n", i, startingPrimeIndex, endingPrimeIndex);
    }
}

void showIndexes() {
//    std::vector<std::size_t> primesIndexes = {
//            4495744, 26272001, 3661285, 33914092,
//            2856137, 96593623, 750612
//    };
//
//    std::size_t primesTotal = 100'000'000;
//    std::size_t threadsTotal = 256 * 256;
//    std::size_t primesPerThreadTotal = primesTotal / threadsTotal;
//    std::size_t primesPerIteration = 64;
//
//    for(const auto index: primesIndexes) {
//        for(std::size_t threadId = 0; threadId < 256 * 256; ++threadId) {
//            const auto threadStart = threadId * primesPerThreadTotal;
//
//            for(std::size_t threadIteration = 0; threadIteration * primesPerIteration < primesPerThreadTotal; ++threadIteration) {
//                const std::size_t iterationStart = threadStart + threadIteration * primesPerIteration;
//
//                if(index > iterationStart && index < iterationStart + primesPerIteration) {
//                    std::cout << "Index " << index << " is " << index - iterationStart << "'s in thread " << threadId << ", iteration " << threadIteration << std::endl;
//                    break;
//                }
//            }
//        }
//    }
}

//void makeTest() {
//    //4495640, ending prime: 4495704
//    const auto primes = loadPrimes("../../all-primes-32-bit.bin");
//
//    Aeu<4096> base = Aeu<4096>(2) * 3 * 7 * 83 * 257;
//
//    for (std::size_t i = 4495704; i < 4495704 + 64; ++i) {
//        base *= primes[i];
//    }
////    base *= primes[4495744];
//    std::cout << "Base: " << base << std::endl;
//
//    // 2 × 3 × 7 × 83 × 257 × 76840783
//
//    const Aeu<4096> a = 2, n = "0x66ee64b74c3d88a79c57064d8365";
//
//    const auto powm = Aeu<4096>::powm(a, base, n);
//    const auto gcd = Aeu<4096>::gcd(powm - 1, n);
//
//    std::cout << "Powm: " << powm << std::endl;
//    std::cout << "GCD: " << gcd << std::endl;
//}




Uns countE(unsigned B, const unsigned* const primes, std::size_t primesCount) {
    auto primeUl = primes[0];

    Uns e = 1;
    for (unsigned pi = 0; primeUl < B && pi < primesCount; ++pi) {
        const auto power = static_cast<unsigned>(log(static_cast<double>(B)) / log(static_cast<double>(primeUl)));
        const auto factor = static_cast<unsigned>(pow(static_cast<double>(primeUl), static_cast<double>(power)));

        std::cout << " * " << primeUl << '^' << power;
//        e *= factor;
        primeUl = primes[pi + 1];
    }
    std::cout << std::endl;

    return e;
}

std::vector<Uns> loadPowers(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load powers table: bad input file");

    const auto recordsNumber = static_cast<std::size_t>(std::filesystem::file_size(fromLocation) / (static_cast<std::size_t>(Uns::getBitness() / 8)));
    std::vector<Uns> records (recordsNumber);

    std::ifstream input(fromLocation, std::ios::binary);
    for(auto& value: records)
        value.readBinary(input);

    return records;
}


template <typename Elem>
int binary_search_find_index(std::vector<Elem> v, Elem data) {
    auto it = std::lower_bound(v.begin(), v.end(), data);
    if (it == v.end()) {
        return -1;
    } else {
        return std::distance(v.begin(), it);
    }
}

using primeType = unsigned;
std::vector<primeType> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");

    const auto primesCount = std::filesystem::file_size(fromLocation) / sizeof(primeType);
    std::ifstream input(fromLocation, std::ios::binary);

    std::vector<primeType> primes (primesCount);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(primeType));

    return primes;
}



std::vector<primeType> loadPrimesFilterBitness(const std::filesystem::path& primesLocation, std::size_t requiredBitness) {
    constexpr auto getFirstWithBitness = [] (std::size_t bitness) {
        if(bitness > sizeof(primeType) * 8)
            throw std::invalid_argument("Specified bitness overflows the type");
        return static_cast<primeType>(1ull << (bitness - 1));
    };
    constexpr auto getLastWithBitness = [] (std::size_t bitness) {
        if(bitness > sizeof(primeType) * 8)
            throw std::invalid_argument("Specified bitness overflows the type");
        return static_cast<primeType>((1ull << bitness) - 1);
    };

    const auto allPrimes = loadPrimes(primesLocation);
    const auto first = getFirstWithBitness(requiredBitness),
        last = getLastWithBitness(requiredBitness);
    const auto firstIdx = binary_search_find_index(allPrimes, first),
        lastIdx = binary_search_find_index(allPrimes, last);

    return { allPrimes.begin() + firstIdx - 1, allPrimes.begin() + lastIdx + 1};
}

constexpr auto getNumbersBitness = [] (primeType number) {
    return sizeof(primeType) * 8 - std::countl_zero(number);
};

int main() {
//    makeTest();
//    showIndexes();
//    kernel(100'000'000, 128);
//         203'227'136

//    //0x104e822c30254e465a6cb92d99617bbca06b53f5d200
//    using Uns = Aeu<2272>;
//    Uns value = "0x104e822c30254e465a6cb92d99617bbca06b53f5d200";

//    std::cout << "{ ";
//    for(auto block: value.blocks)
//        std::cout << "0x" << std::hex << block << ", ";

//    std::array<unsigned, 71> blocks = {0x53f5d200, 0x7bbca06b, 0xb92d9961, 0x4e465a6c, 0x822c3025, 0x104e, 0x0, 0x0,
//                                        0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
//                                        0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
//                                        0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
//                                        0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
//    };
//
//    Uns value2;
//    for(std::size_t i = 0; i < value2.blocks.size(); ++i)
//        value2.blocks[i] = blocks[i];
//
//    std::cout << value << std::endl << value2 << std::endl;
//
//    std::cout << value2.blocks.size() << std::endl;

//    Uns base = Uns(1) * 512 * 729 * 3125 * 2401 * 121 *
//               169 * 17 * 19 * 21 * 23 * 27 * 29
//               * 31 * 37 * 83 * 131 * 199 * 257
//               * 317 * 1051 * 1297 * 1637 * 1933;
//
//    std::cout << base << std::endl << base.bitCount() << std::endl;

//    Uns base = Uns(1) * 512 * 729 * 3125 * 2401 * 121 *
//    169 * 17 * 19 * 21 * 23 * 27 * 29
//    * 31 * 37 * 79 * 79 * 83 * 131 * 137 * 199 * 257
//    * 317 * 359 * 601 * 1051 * 1297 * 1637 * 1933 * 24671;
//
//    std::cout << base << std::endl << base.bitCount() << std::endl;

//    const auto powers = loadPowers("../../powers/powers-8192.txt");
//    std::cout << std::hex << "0x" << powers.back() << std::endl;

//    std::cout << primes[3745306] << std::endl;
//    std::cout << binary_search_find_index(primes, 63287041u) << std::endl;

//    std::cout << "Your bitness: " << std::endl;
//    unsigned bitness; std::cin >> bitness;
//    std::cout << "Received bitness " << bitness << std::endl;

//    const std::filesystem::path primesLocation = "../../all-primes-32-bit.bin";
//    const auto values = loadPrimesFilterBitness(primesLocation, 25);
//
//    std::cout << "First: " << values[0] << " = " << std::bitset<32>(values[0]) << " (" << getNumbersBitness(values[0]) << " bits)" << std::endl;
//    std::cout << "Second: " << values[1] << " = " << std::bitset<32>(values[1]) << " (" << getNumbersBitness(values[1]) << " bits)" << std::endl;
//
//    std::cout << "Before last: " << values[values.size() - 2] << " = " << std::bitset<32>(values[values.size() - 2]) << " (" << getNumbersBitness(values[values.size() - 2]) << " bits)" << std::endl;
//    std::cout << "Last: " << values.back() << " = " << std::bitset<32>(values.back()) << " (" << getNumbersBitness(values.back()) << " bits)" << std::endl;

//    Uns base = Uns(1) * 512 * 729 * 3125 * 2401 * 121 *
//               169 * 17 * 19 * 21 * 23 * 27 * 29
//               * 31 * 37 * 83 * 131 * 199 * 257
//               * 317 * 1051 * 1297 * 1637 * 1933;
//    std::cout << std::hex << "0x" << base << std::endl;

//    std::cout << getNumbersBitness(63287041) << std::endl;

//    Uns base = Uns(2) * 3 * 5 * 17 * 19 * 21 * 23 * 27 * 29
//               * 31 * 37 * 79 * 79 * 83 * 121 * 131 * 137
//               * 169 * 199 * 257 * 317 * 359 * 512 * 601 * 729 * 1051 * 1297 * 1637 * 1933 * 2401 * 3125 * 24671;
//    std::cout << std::hex << "0x" << base << std::endl;


}