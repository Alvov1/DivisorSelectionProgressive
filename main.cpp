#include <iostream>
#include <fstream>
#include <filesystem>
#include <primesieve.hpp>
#include <thrust/device_vector.h>
#include "AesiMultiprecision.h"

std::vector<uint64_t> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::exists(fromLocation) || !std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");

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
        throw std::runtime_error("Failed to save prime table: bad output file");

    const uint64_t primesCount = primes.size();
    output.write(reinterpret_cast<const char*>(&primesCount), sizeof(uint64_t));

    for(const auto& prime: primes)
        output.write(reinterpret_cast<const char*>(&prime), sizeof(uint64_t));
}

thrust::device_vector<uint64_t> deviceLoadPrimes(const std::filesystem::path& fromLocation) {
    return loadPrimes(fromLocation);
}

//gpu void kernel(const Aesi& n, const thrust::device_vector<uint64_t>& primes) {
//    const auto threadIdx = blockDim.x * blockIdx.x + threadIdx.x,
//        threads = gridDim.x * blockDim.x,
//        max_it = 400000 / threads,
//        bStart = 2 + blockIdx.x,
//        bInc = gridDim.x,
//        B_MAX = 2000000000;
//
//    const auto checkWriteRepeat = [&n] (const Aesi<512>& value) {
//        if(value < 2 || value >= n) return false;
//        char buffer[512] {}; value.getString<10>(buffer, 512);
//        return printf("Found divisor: %s\n", buffer) > 14;
//    };
//
//    Aesi a = threadIdx * max_it + 2, e = 1;
//    for(unsigned B = bStart; B < B_MAX; B += bInc) {
//        const auto primeUl = primes[0];
//
//        for(unsigned pi = 0; primeUl < B; ++pi) {
//            const unsigned power = log(static_cast<double>(B)) / log(static_cast<double>(primeUl));
//            e *= static_cast<uint64_t>(pow(static_cast<double>(primeUl), static_cast<double>(power)));
//            primeUl = primes[pi + 1];
//        }
//
//        if(e == 1) continue;
//
//        for(unsigned it = 0; it < max_it; ++it) {
//            if(checkWriteRepeat(Aesi<512>::gcd(a, n)))
//                return;
//
//            if(checkWriteRepeat(Aesi<512>::gcd(Aesi<512>::powm(a, e, n) - 1, n)))
//                return;
//
//            a += threads * max_it;
//        }
//    }
//}

gpu void kernel(const Aesi& value) {
    const auto threadIdx = blockDim.x * blockIdx.x + threadIdx.x;
    if(threadIdx > 0) return;

    char buffer[512] {}; value.getString<10>(buffer, 512);
    printf("Kernel thread: factorizing value %s\n", buffer);
}

int main(int argc, const char* const* const argv) {
    using namespace std::string_view_literals;
    if(argc < 3)
        return std::printf("Usage:\n\t%s load-primes <primes location>"
                           "\n\t%s generate-primes <primes count> <primes location>"
                           "\n\t%s factorize <number> <primes location>", argv[0], argv[0], argv[0]);

    try {
        if (argc == 3 && argv[1] == "load-primes"sv) {
            const auto primes = loadPrimes(argv[2]);
            std::cout << "Loaded prime table of " << primes.size() << " elements." << std::endl;
        } else if (argc > 3) {
            if (argv[1] == "generate-primes"sv) {
                std::vector<uint64_t> primes(std::stoi(argv[2]));
                primesieve::iterator it;
                for (auto &prime: primes)
                    prime = it.next_prime();
                savePrimes(primes, argv[3]);
                std::cout << "Generated prime table of " << primes.size() << " elements to '" << argv[3] << "'." << std::endl;
            } else if (argc > 3 && argv[1] == "factorize"sv) {
                Aesi number = std::string_view(argv[2]);
                std::cout << "Factorizing number " << std::hex << std::showbase << number << '.' << std::endl;

//                const auto primeTable = deviceLoadPrimes(argv[3]);
//                std::cout << "Primes are loaded to device" << std::endl;

                kernel<<<32, 32>>>(number);
                if (cudaSuccess != cudaDeviceSynchronize())
                    throw std::runtime_error("Kernel launch failed");
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Failed: " << e.what() << '.' << std::endl;
    }

    return 0;
}
