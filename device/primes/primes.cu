#include <iostream>
#include <fstream>
#include <filesystem>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "Aeu.h"
#include "Timer.h"
#include "primes-kernel.h"

using primeType = unsigned;
thrust::host_vector<primeType> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");

    const auto primesCount = std::filesystem::file_size(fromLocation) / sizeof(primeType);
    std::ifstream input(fromLocation, std::ios::binary);

    std::vector<primeType> primes (primesCount);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(primeType));

    return primes;
}


int main(int argc, const char* const* const argv) {
    if(argc < 5)
        return std::printf("Usage: %s <number> <primes location> <threads> <iterations>", argv[0]);

    const Uns number = std::string_view(argv[1]);
    thrust::device_vector<Uns> numberAndFactor = { number, Uns { 0 } };
    Timer::init() << "Factorizing number " << std::hex << std::showbase << number << std::dec << " (" << number.bitCount() << " bits)." << Timer::endl;

    const thrust::device_vector<primeType> primes = loadPrimes(argv[2]);
    Timer::out << "Loaded table of primes with " << primes.size() << " elements." << Timer::endl;

    const auto threads = std::stoul(argv[3]), iterations = std::stoul(argv[4]);
    const auto timePoint = std::chrono::system_clock::now();
    const auto timeT = std::chrono::system_clock::to_time_t(timePoint);
    Timer::out << std::ctime(&timeT) << " Starting kernel <<<" << threads << ", " << threads << ">>>. Using bitness " << Uns::getBitness() << ". Iterations: " << iterations << Timer::endl;
    kernel<<<threads, threads>>>(
            thrust::raw_pointer_cast(numberAndFactor.data()),
            thrust::raw_pointer_cast(primes.data()),
            primes.size(),
            iterations);

    const auto code = cudaDeviceSynchronize();
    if(code != cudaSuccess)
        return std::printf("Kernel launch failed: %s.\n", cudaGetErrorString(code));

    const Uns factor = numberAndFactor[1];
    if(factor != 0 && number % factor == 0)
        Timer::out << "Kernel completed. Founded factor: " << std::hex << std::showbase << factor << '.' << Timer::endl;
    else Timer::out << "Kernel completed, but factor " << std::hex << std::showbase << factor << " is incorrect." << Timer::endl;

    return 0;
}