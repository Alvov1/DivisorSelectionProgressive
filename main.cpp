#include <iostream>
#include <fstream>
#include <filesystem>
#include <primesieve.hpp>
#include "AesiMultiprecision.h"

gpu void kernel() {

}

template <typename Integral = uint64_t> requires (std::is_integral_v<Integral>)
std::vector<Integral> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::exists(fromLocation) || !std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file.");

    std::ifstream input(fromLocation, std::ios::binary);
    Integral buffer {}; input.read(reinterpret_cast<char*>(&buffer), sizeof(Integral));

    std::vector<Integral> primes (buffer);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(Integral));

    return primes;
}

template <typename Integral = uint64_t> requires (std::is_integral_v<Integral>)
void savePrimes(const std::vector<Integral>& primes, const std::filesystem::path& toLocation) {
    std::ofstream output(toLocation, std::ios::binary);
    if(output.fail())
        throw std::runtime_error("Failed to save prime table: bad output file.");

    const uint64_t primesCount = primes.size();
    output.write(reinterpret_cast<const char*>(&primesCount), sizeof(Integral));

    for(const auto& prime: primes)
        output.write(reinterpret_cast<const char*>(&prime), sizeof(Integral));
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
    }

    return 0;
}
