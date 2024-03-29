#include <iostream>
#include <fstream>
#include <filesystem>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "Aesi.h"
#include "Timer.h"
#include "primes-kernel.h"

thrust::host_vector<uint64_t> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");

    std::ifstream input(fromLocation, std::ios::binary);
    uint64_t buffer {}; input.read(reinterpret_cast<char*>(&buffer), sizeof(uint64_t));

    thrust::host_vector<uint64_t> primes (buffer);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(uint64_t));

    return primes;
}


int main(int argc, const char* const* const argv) {
    if(argc < 5)
        return std::printf("Usage: %s <number> <primes location> <<64> threads parameters <64>>", argv[0]);

    const Uns number = std::string_view(argv[1]);
    thrust::device_vector<Uns> numberAndFactor = { number, Uns { 0 } };
    Timer::init() << "Factorizing number " << std::hex << std::showbase << number << std::dec << " (" << number.bitCount() << " bits)." << Timer::endl;

    const thrust::device_vector<uint64_t> primes = loadPrimes(argv[2]);
    Timer::out << "Loaded table of primes with " << primes.size() << " elements." << Timer::endl;

    const auto lThreads = std::stoi(argv[3]), rThreads = std::stoi(argv[4]);
    Timer::out << "Starting kernel <<<" << lThreads << ", " << rThreads << ">>>. Using bitness " << Uns::getBitness() << '.' << Timer::endl;
    kernel<<<lThreads, rThreads>>>(
            thrust::raw_pointer_cast(numberAndFactor.data()),
            thrust::raw_pointer_cast(primes.data()),
            primes.size());

    const auto code = cudaDeviceSynchronize();
    if(code != cudaSuccess)
        return std::printf("Kernel launch failed: %s.\n", cudaGetErrorString(code));

    const Uns factor = numberAndFactor[1];
    if(factor != 0 && number % factor == 0)
        Timer::out << "Kernel completed. Founded factor: " << std::hex << std::showbase << factor << '.' << Timer::endl;
    else Timer::out << "Kernel completed, but factor " << std::hex << std::showbase << factor << " is incorrect." << Timer::endl;

    return 0;
}