#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <bitset>
#include <Aeu.h>


using Uns = Aeu<9216>;

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


int binary_search_find_index(const std::vector<unsigned>& v, unsigned data) {
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

    //30'000;
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
    const std::filesystem::path primesLocation = "../../all-primes-32-bit.bin";
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

    //2^2 × 3^3 × 5 × 97 × 23747 × 105'137 × 745'335'889
    //2^2 × 3^3 × 5 × 97 × 23747 × 105'137 × 745'335'889
//    const auto values = loadPrimes(primesLocation);
//    Uns value = Uns(18) * 23747 * 105'137;// * 105'137;
////    Uns value = Uns(2) * 2 * 3 * 3 * 3 * 5 * 97 * 23747 * 105'137 * 745'335'889;
//    for(auto tValue: values) {
//        value *= tValue;
//        if(tValue > 3000)
//            break;
//    }
//
//    value *= 745'335'889;
    //2 * 3 × 5 × 29 × 79^2 × 137 × 359 × 601 × 24671 × 63287041
//    std::cout << std::hex << "0x" << value << " (" << std::dec << value.bitCount() << ")" << std::endl;


//    Aeu<2048> v = Aeu<2048>(1) * 2 * 3 * 5 * 29 * 79 * 79 * 137 * 359 * 601 * 24671 * 63287041;
//
//
    const auto primes = loadPrimesFilterBitness(primesLocation, 30);
    std::cout << primes.size() << std::endl;
    std::cout << binary_search_find_index(primes, 745'335'889) << std::endl;

//    Uns n = "0x2c78e5029bd4bb31f69c965cb016b622c0540e85e919ed563f27a21098084950addf";
//    Aeu<4096> v = "0x7b44126612c8e50fc5fed536db8b3f4e955b023bfe2a8f55b9c7ef82d47fbad084cff562269114488d9596cae7277ac7113560d4b9ce68b3142e5789884ee1b9e8aad69359bffb4258827525a03ffca96c59be362dce6c7f107db867f9d268d052bc42b9197496e05c471bfe2303a2287923d4e13713518e242b04d9597e911d9d9ba302d79193c63d48c510ee9ffc84495b8a1797b0616a4d69a09ca16e2dafa34b0250473b1389dba16479c2bd2a4b46fa01d9cef5c321dab59bcab53ed2439a6ec708f447e80790392a71e3b13df6cd1dccb2b03bd574a1b56919c524202f84fd5b2bf4563e50cba9cde095ab376b8ca344954a8f11706530adf42a5b0033061fbf9c86b9bfedefbaa160aa6e17c851bede26ece39fcde0c0eb213266868559c4c93624bb8351d4889223a8f6b02fa0c8122c52e5c96e1cf22c0a78d12dda1754dac66cb954c8ab248d99e14600f3e89fb9b059dac625583164b6ccbb95391eac9f2944b8db27482cf14bf42686ee4db5d0e3745fa15f0b04e58789184e1ac66d1eac533786547c665d1a32ea7c14f3ed31047c4003fee5d72036550c6c396633b885cebf95107254073bb06f462b4b49349ea69c73705140735bb1b213d8fdbb34f8aa8e31ff25f1734cc2cfb92301414bce2ce122b0c47cbb39ff49e58ac2b30b8d151af8cfa8fec7c6f034c9a10e146f43df610ab6767893b07509193ee15bbacd1028f2c2997d69d6c78bce1d34e995e88a";
//    value *= primes[1681618];
//    value *= 63287041;



//    std::cout << Uns::gcd(Uns::powm(2, value, n) - 1, n) << std::endl;

//    1894122 -> 1681618
//    1681596
//    1681618
//    const Aeu<4096> a = 2, n = "0x545e18a29cdbb180605b5ab2e8a7012948e49357";
//
//    Aeu<4096> v = "0x79df1b14b1ea3d0f65603047b7e2627c007e944f01e6e1d02d28d612530cf679407ac199ad956e454fbcd42285c1780f2173265831de519ae090ec76a0cbd93cb1cbc2cb79a477b969b11e527efea891ca9492e549e08999b5567f1f380d27d1e92ea80b5a41c52316565c43fa5298ebc792c31eee7136e12f56d0c90eb08f7da484691550a9d8db9e7bd8a109e1281616dcd59fd79f0a6aad77401f82a7af6e1094782a5baba0f49a881cb1003695664933c733150176682b620278a2b2b040d494c92987264b33132c7a6faea20dfb645409d8b57e966d9bc69285225379506e58ed3b7f713ff282b30d48c8267a0f96b4066312b4587eb6a2067f06fe84a10321d5e6e05a8a593089ec0eb29cdc8ae78b2183b84d6526e05cfbec66292f559169a23be";
////    for (unsigned i = 0; i < 28; ++i) {
////        v *= 63287041;
////    }
//
//    v *= 652450049;
//    v *= 2;
//
//    std::cout << v.bitCount() << ". Multiplied. " << v.bitCount() << std::endl;
//
//    std::cout << Aeu<4096>::gcd(Aeu<4096>::powm(a, v, n) - 1, n);
//
////    Aeu<2048> v = Aeu<2048>(1) * 2 * 3 * 5 * 29 * 79 * 79 * 137 * 359 * 601 * 24671 * 63287041;
////        std::cout << Aeu<2048>::gcd(Aeu<2048>::powm(a, v, n) - 1, n);


    /* Итерация по элементам, с шагом по количеству потоков. */
    for(std::size_t i = threadId; i < primesCount; i += threadsCount) {
        /* В стеке хранится множитель, и его позиция в общем массиве. */
        PrimeStack<primeType, stackSize + 2> stack {};
        std::size_t primesPosition = i;

        do {
            /* Если стек уже заполнен: проверка текущего произведения. */
            if (stack.getSize() >= stackSize) {
                const Uns product = stack.unite<Uns>() * startingBase * primes[i];
                const auto candidate = Uns::gcd(Uns::powm(a, product, n) - 1, n);
                if(gcd > 1u)
                    return factor.tryAtomicSet(gcd);
            }

            /* Извлечение элементов из стека. */
            if (stack.getSize() >= stackSize || primesPosition >= primesCount) {
                std::size_t checkedLast = primesCount;
                do {
                    primesPosition = stack.pop();
                    checkedLast -= primesPositionShift;
                } while (primesPosition == checkedLast);

                primesPosition += primesPositionShift;
                stack.push(primes[primesPosition], primesPosition);
                primesPosition += primesPositionShift;
            } else {
                /* Добавление новых элементов. */
                stack.push(primes[primesPosition], primesPosition);
                primesPosition += primesPositionShift;
            }

        } while (!stack.isEmpty());
    }
}