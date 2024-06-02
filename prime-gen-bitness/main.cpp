#include <iostream>
#include <fstream>
#include <primesieve.hpp>
#include <limits>
#include <vector>
#include <filesystem>

using primeType = uint64_t;
std::size_t binary_search_find_index(const std::vector<primeType>& v, primeType data) {
    auto it = std::lower_bound(v.begin(), v.end(), data);
    if (it == v.end()) {
        return 0;
    } else {
        return std::distance(v.begin(), it);
    }
}


std::vector<primeType> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");


    const auto primesCount = std::filesystem::file_size(fromLocation) / sizeof(primeType);
    std::cout << "Loading prime table of " << primesCount << " elements ("
              << std::filesystem::file_size(fromLocation) << " bytes)" << std::endl;
    std::ifstream input(fromLocation, std::ios::binary);

    std::vector<primeType> primes (primesCount);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(primeType));

    return primes;
}

template <typename Type>
std::size_t getBitness (Type prime) {
    return sizeof(Type) * 8 - std::countl_zero(prime);
}

int main(int argc, const char* const* const argv) {
    const auto primes = loadPrimes("../../primes-34bit-only.bin");
    const auto value = 14'786'747'089ULL;
    std::cout << "Value: " << value << ", index: " << binary_search_find_index(primes, value);
    return 0;
}