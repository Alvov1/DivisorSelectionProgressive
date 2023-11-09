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

    if(argc < 2 || argv[1] == "generate-primes"sv) {
        std::vector<uint64_t> primeTable(80'000'000); //78'643'200
        primesieve::iterator it;
        for(unsigned i = 0; i < primeTable.size(); ++i) {
            primeTable[i] = it.next_prime();
            if((i + 1) % 1'000'000 == 0)
                std::cout << "Completed " << i + 1 << "." << std::endl;
        }

        savePrimes(primeTable, "primes.txt");
    } else if(argv[1] == "load-primes"sv) {
        const auto primes = loadPrimes("primes.txt");
        std::cout << "Loaded prime table of " << primes.size() << " elements." << std::endl;
    }

    return 0;
}
