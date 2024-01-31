#include <iostream>
#include <fstream>
#include <filesystem>
#include "Aesi.h"

std::vector<uint64_t> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");

    std::ifstream input(fromLocation, std::ios::binary);
    uint64_t buffer {}; input.read(reinterpret_cast<char*>(&buffer), sizeof(uint64_t));

    std::vector<uint64_t> primes (buffer);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(uint64_t));

    return primes;
}

template<std::size_t precision>
Aesi<precision> countE(unsigned B, const std::vector<uint64_t>& primes) {
    auto primeUl = primes[0];

    Aesi<precision> e = 1;

    for (unsigned pi = 0; primeUl < B && pi < primes.size(); ++pi) {
        const auto power = static_cast<unsigned>(log(static_cast<double>(B)) / log(static_cast<double>(primeUl)));
        e *= static_cast<uint64_t>(pow(static_cast<double>(primeUl), static_cast<double>(power)));
        primeUl = primes[pi + 1];
    }

    return e;
}

template <std::size_t precision>
void preparePowers(const std::vector<uint64_t>& primes, const std::filesystem::path& outputLocation) {
    Aesi<precision> previous = 0;

    std::ofstream output(outputLocation, std::ios::binary);

    for (unsigned B = 2; ; B++) {
        const auto e = countE<precision>(B, primes);

        if(e > previous) {
            if(e.bitCount() >= 32) {
                std::cout << "B: " << std::dec << B << ", e: 0x" << std::hex << e << ", " << std::dec << e.bitCount() << " bits" << std::endl;
                output.write(reinterpret_cast<const char*>(&B), sizeof(unsigned));
                e.writeBinary(output);
            }

            previous = e;

            if(e.bitCount() >= static_cast<std::size_t>(Aesi<precision>::getBitness() * 0.99))
                break;
        }
    }
}

int main() {
    const std::filesystem::path primes = "../../primes.txt", output = "../../powers/powers-1024-tuples.txt";
    preparePowers<1024>(loadPrimes(primes), output);

    return 0;
}