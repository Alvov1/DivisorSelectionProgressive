#include <iostream>
#include <fstream>
#include <primesieve.hpp>
#include <limits>

int main(int argc, const char* const* const argv) {
    using namespace std::string_literals;
    if(argc < 3)
        return std::printf("Usage: %s <smooth location> <options>\nOptions: boarder 32 bit number | full\n", argv[0]);

    try {
        std::ofstream output(argv[1], std::ios::binary);
        if (output.fail())
            return std::printf("Failed to create the output file %s\n", argv[1]);

        if(argv[2] == "full"s) {
            primesieve::iterator it;
            for (uint64_t prime = it.next_prime();
                 prime < std::numeric_limits<unsigned>::max(); prime = it.next_prime()) {
                const auto casted = static_cast<unsigned>(prime);
                output.write(reinterpret_cast<const char *>(&casted), sizeof(unsigned));
            }
        } else {
            const auto boarder = std::stoull(argv[2]);
            primesieve::iterator it;
            for (uint64_t prime = it.next_prime(); prime < boarder; prime = it.next_prime()) {
                const auto casted = static_cast<unsigned>(prime);
                output.write(reinterpret_cast<const char *>(&casted), sizeof(unsigned));
            }
        }

        std::cout << "Generated prime table '" << argv[1] << "'." << std::endl;
    } catch (const std::exception& e) {
        return std::printf("Failed: %s\n", e.what());
    }
    return 0;
}