#include <iostream>
#include <fstream>
#include <filesystem>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "Aeu.h"
#include "Timer.h"

using primeType = unsigned;
thrust::host_vector<primeType> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");

    const auto primesCount = std::filesystem::file_size(fromLocation) / (sizeof(primeType) * 2);
    std::ifstream input(fromLocation, std::ios::binary);

    std::vector<primeType> primes (primesCount);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(primeType));

    return primes;
}








using ValuePrecision = Aeu128;
using TestPrecision = Aeu<ValuePrecision::getBitness() + 32 * 4>;

__global__
void kernel(ValuePrecision* const numberAndFactor,
            const primeType* const primes,
            std::size_t primesCount) {
    const unsigned threadId = blockDim.x * blockIdx.x + threadIdx.x;
    if(threadId % 32768 != 0) return;

    char buffer [512] {};
    e.getString<16>(buffer, 512, true, false);

    printf("Thread %u. 0x%s\n", threadId, buffer);
}











int main(int argc, const char* const* const argv) {
    try {
        if(argc < 5)
            return std::printf("Usage: %s <number> <smooth location> <threads> <iterations>", argv[0]);

        const ValuePrecision number = std::string_view(argv[1]);
        thrust::device_vector<ValuePrecision> numberAndFactor = { number, ValuePrecision { 0 } };
        Timer::init() << "Factorizing number 0x" << std::hex << number << std::dec << " (" << number.bitCount() << " bits)." << Timer::endl;

        const thrust::device_vector<primeType> primes = loadPrimes(argv[2]);
        Timer::out << "Loaded table of primes with " << primes.size() << " elements." << Timer::endl;

        const auto threads = std::stoul(argv[3]), iterations = std::stoul(argv[4]);
        const auto timePoint = std::chrono::system_clock::now();
        const auto timeT = std::chrono::system_clock::to_time_t(timePoint);
        Timer::out << std::ctime(&timeT) << " Starting kernel <<<" << threads << ", " << threads << ">>>. Using bitness "
            << TestPrecision::getBitness() << ". ITERATIONS PER PRIME: " << iterations << Timer::endl;
        kernel<<<threads, threads>>>(
                thrust::raw_pointer_cast(numberAndFactor.data()),
                thrust::raw_pointer_cast(primes.data()),
                primes.size());

        const auto code = cudaDeviceSynchronize();
        if(code != cudaSuccess)
            return std::printf("Kernel launch failed: %s.\n", cudaGetErrorString(code));

        const ValuePrecision factor = numberAndFactor[1];
        if(factor != 0 && number % factor == 0)
            Timer::out << "Kernel completed. Founded factor: 0x" << std::hex << factor << '.' << Timer::endl;
        else Timer::out << "Kernel completed, but factor 0x" << std::hex << factor << " is incorrect." << Timer::endl;
    } catch (const std::exception& e) {
        return std::printf("Failed: %s.\n", e.what());
    }
}