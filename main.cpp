#include <iostream>
#include <fstream>
#include <filesystem>
#include <primesieve.hpp>
#include "AesiMultiprecision.h"

std::vector<uint64_t> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::exists(fromLocation) || !std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file.");

    std::ifstream input(fromLocation, std::ios::binary);
    uint64_t buffer {}; input.read(reinterpret_cast<char*>(&buffer), sizeof(uint64_t));

    std::vector<uint64_t> primes (buffer);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(uint64_t));

    return primes;
}

void savePrimes(const std::vector<uint64_t>& primes, const std::filesystem::path& toLocation) {
    std::ofstream output(toLocation, std::ios::binary);
    if(output.fail())
        throw std::runtime_error("Failed to save prime table: bad output file.");

    const uint64_t primesCount = primes.size();
    output.write(reinterpret_cast<const char*>(&primesCount), sizeof(uint64_t));

    for(const auto& prime: primes)
        output.write(reinterpret_cast<const char*>(&prime), sizeof(uint64_t));
}

gpu void kernel() {

}


int main(int argc, const char* const* const argv) {
    using namespace std::string_view_literals;
    if(argc < 3)
        return std::printf("Usage:\n\t%s generate-primes <primes location>\n\t%s load-primes <primes location>", argv[0], argv[0]);

    if(argv[1] == "generate-primes"sv) {
        std::vector<uint64_t> primes(134'217'727);
        primesieve::iterator it;
        for(auto& prime: primes)
            prime = it.next_prime(); savePrimes(primes, argv[2]);
        std::cout << "Generated prime table of " << primes.size() << " elements." << std::endl;
    } else if(argv[1] == "load-primes"sv) {
        const auto primes = loadPrimes(argv[2]);
        std::cout << "Loaded prime table of " << primes.size() << " elements." << std::endl;
    } else if(argv[1] == "factorize"sv) {
        Aesi number = std::string_view(argv[2]);
        std::cout << "Factorizing number " << std::hex << number << '.' << std::endl;
    }

    return 0;
}
