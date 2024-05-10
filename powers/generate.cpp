#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <unordered_map>
#include <map>
#include "Aesi.h"

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

template<std::size_t precision>
Aeu<precision> countE(unsigned B, unsigned shift, const std::vector<uint64_t>& primes) {
    auto primeUl = primes[0];

    Aeu<precision> e = 1;

    std::cout << "B: " << B << ", value: ";

    for (unsigned pi = 0; primeUl < B && pi < primes.size(); ++pi) {
        auto power = static_cast<unsigned>(log(static_cast<double>(B)) / log(static_cast<double>(primeUl)));
        if(power > shift)
            power -= shift;

        std::cout << primeUl << '^' << power << " * ";
        e *= { static_cast<uint64_t>(pow(static_cast<double>(primeUl), static_cast<double>(power))) };
        primeUl = primes[pi + 1];
    }

    std::cout << std::endl;

    return e;
}

template <std::size_t precision>
void preparePowers(const std::vector<uint64_t>& primes, std::ostream& output, unsigned shift) {
    Aeu<precision> previous = 0;

    unsigned B = 0;
    for (B = 2; ; B += 25) {
        for (unsigned i = 0; i < 5; ++i) {
            const auto e = countE<precision>(B, i, primes);

            if(e > previous) {
                /*if(e.bitCount() >= 32) {
                    std::cout << "B: " << std::dec << B << ", e: 0x" << std::hex << e << ", " << std::dec << e.bitCount() << " bits" << std::endl;
                    output.write(reinterpret_cast<const char*>(&B), sizeof(unsigned));
                    e.writeBinary(output);
                } */

                previous = e;

                if(e.bitCount() >= static_cast<std::size_t>(Aesi<precision>::getBitness() * 0.99)) {
                    std::cout << std::endl << "Last: " << previous << std::endl;
                    return;
                }
            }
        }

        std::cout << std::endl;
    }


}

template <std::size_t precision>
void preparePowers(const std::vector<uint64_t>& primes, std::ostream& output) {

    /* 1. At first, we need to find the maximum B for which product overflows the precision. */
    auto [factors, B] = [&primes] {
        std::pair values = { std::map<uint64_t, uint8_t>{}, unsigned() };
        auto& [factorsWithPowers, inB] = values;

        for(inB = 48; inB += 5;) {
            Aeu<precision + 128> buffer = 1u;

            for(const auto& prime: primes) {
                if(prime >= inB)
                    break;

                const auto power = static_cast<uint8_t>(log(static_cast<double>(inB)) / log(static_cast<double>(prime)));
                buffer *= { static_cast<uint64_t>(pow(static_cast<double>(prime), static_cast<double>(power))) };
                if(buffer.bitCount() > precision)
                    return values;

                factorsWithPowers[prime] = power;
            }
        }

        return values;
    } ();
    [&factors, &B] {
        Aeu<precision> buffer = 1;
        for(auto& [factor, power]: factors) {
            std::cout << factor << '^' << static_cast<unsigned>(power) << " * ";
            buffer *= { static_cast<uint64_t>(pow(factor, power)) };
        }
        std::cout << std::endl << buffer << ". Bitness: " << buffer.bitCount() << ". B = " << B << std::endl << std::endl;
    } ();

}



int main() {
    const std::filesystem::path primes = "../../all-primes-32-bit.bin", outputLocation = "../test-test.txt";
//    if(std::filesystem::is_regular_file(outputLocation))
//        return std::printf("Output file exists and no overwritting.");

//    std::ofstream output(outputLocation, std::ios::binary);
//    if(output.fail())
//        return std::printf("Failed to open output file.");

//    preparePowers<96>(loadPrimes(primes), output);

    auto primesL = loadPrimes(primes);
    std::cout << primesL.back() << std::endl;

    return 0;
}