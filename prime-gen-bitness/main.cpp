#include <iostream>
#include <fstream>
#include <primesieve.hpp>
#include <limits>

template <typename Type>
std::size_t getBitness (Type prime) {
    return sizeof(Type) * 8 - std::countl_zero(prime);
}

int main(int argc, const char* const* const argv) {
    using namespace std::string_literals;
//    if(argc < 3)
//        return std::printf("Usage: %s <smooth location> <options>\nOptions: boarder 32 bit number | full\n", argv[0]);

    const auto bitness = 34;

    try {
        std::size_t count = 1;

        primesieve::iterator it;
        for (uint64_t prime = it.next_prime();
             getBitness(prime) < bitness; prime = it.next_prime());

        std::cout << "Started." << std::endl;
        for(uint64_t prime = it.next_prime();
            getBitness(prime) < bitness + 1; prime = it.next_prime())
            count += 1;

        std::cout << "Bitness: " << bitness << " - " << count << " primes (" << static_cast<double>(count * 8) / (1024. * 1024. * 1024.) << " GBytes)" << std::endl;
    } catch (const std::exception& e) {
        return std::printf("Failed: %s\n", e.what());
    }
    return 0;
}