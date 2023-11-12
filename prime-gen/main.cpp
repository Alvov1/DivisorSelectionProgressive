#include <iostream>
#include <fstream>
#include <primesieve.hpp>
#include "Timer.h"

int main(int argc, const char* const* const argv) {
    if(argc < 3)
        return std::printf("Usage: %s <primes count> <primes location>\n", argv[0]);

    try {
        std::ofstream output(argv[2], std::ios::binary);
        if (output.fail())
            return std::printf("Failed to create the output file %s\n", argv[2]);

        Timer::init();

        const uint64_t primeCount = std::stoi(argv[1]);
        output.write(reinterpret_cast<const char *>(&primeCount), sizeof(uint64_t));

        primesieve::iterator it;
        for (uint64_t i = 0, prime = 0; i < primeCount; ++i) {
            prime = it.next_prime();
            output.write(reinterpret_cast<const char *>(&prime), sizeof(uint64_t));
        }

        Timer::out << "Generated prime table '" << argv[2] << "' of " << primeCount << " elements." << Timer::endl;
    } catch (const std::exception& e) {
        return std::printf("Failed: %s\n", e.what());
    }

    std::vector<unsigned>().size()
    return 0;
}