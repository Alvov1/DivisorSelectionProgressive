#include <iostream>
#include <fstream>
#include <filesystem>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "Aesi.h"
#include "Timer.h"
#include "powers-kernel.h"

template <std::size_t bitness>
thrust::host_vector<Aesi<bitness>> loadPowers(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load powers table: bad input file");

    const std::size_t fSize = std::filesystem::file_size(fromLocation), oneRecord = bitness / 8, recordsNumber = fSize / (oneRecord);
    if(fSize % oneRecord != 0)
        throw std::invalid_argument("Possible precision misplace.");
    std::ifstream input(fromLocation, std::ios::binary);

    thrust::host_vector<Aesi<bitness>> records (recordsNumber);
    std::size_t previousBitness = 0;
    for(auto& value: records) {
        value.readBinary(input);
        if(previousBitness > value.bitCount())
            throw std::runtime_error("File is corrupted.");
        previousBitness = value.bitCount();
    }

    return records;
}

template<std::size_t bitness>
void launch(std::size_t threads, std::size_t iterations, const char* const argument) {
    using Uns = Aesi<bitness>;

    const Uns number = std::string_view(argument);
    thrust::device_vector<Uns> factorize = { number, Uns { 0 } };
    Timer::out << "Factorizing number " << std::hex << std::showbase << number << std::dec << " (" << number.bitCount() << " bits)." << Timer::endl;

    const auto powersLocation = std::filesystem::current_path().parent_path().parent_path() / "powers" / ("powers-" + std::to_string(bitness) + ".txt");
    const thrust::device_vector<Uns> table = loadPowers<bitness>(powersLocation);
    Timer::out << "Loaded table with " << table.size() << " elements." << Timer::endl;

    kernel<<<threads, threads>>>(
            thrust::raw_pointer_cast(factorize.data()),
            thrust::raw_pointer_cast(table.data()),
            table.size(),
            iterations);

    const auto code = cudaDeviceSynchronize();
    if (code != cudaSuccess)
        throw std::runtime_error("Kernel launch failed: " + std::string(cudaGetErrorString(code)));

    const Uns factor = factorize[1];
    if (factor == 0)
        Timer::out << "Kernel completed, but factor is not found." << Timer::endl;
    else
        if(number % factor == 0)
            Timer::out << "Kernel completed. Founded factor: " << std::hex << std::showbase << factor << '.' << Timer::endl;
        else
            Timer::out << "Kernel completed, but factor is incorrect." << Timer::endl;
}

int main(int argc, const char* const* const argv) {
    try {
        if (argc < 5)
            return std::printf("Usage: %s <number> <bitness> <threads-number> <iterations-count>", argv[0]);

        const std::size_t bitness = std::stoi(argv[2]), threads = std::stoi(argv[3]), iterations = std::stoi(argv[4]);
        Timer::init() << "Using bitness " << bitness << ". Starting threads <<<" << threads << ", " << threads << ">>>. Using iterations " << iterations << Timer::endl;

        if(bitness == 1024)
            launch<1024>(threads, iterations, argv[1]);
        if(bitness == 1536)
            launch<1536>(threads, iterations, argv[1]);
        if(bitness == 2048)
            launch<2048>(threads, iterations, argv[1]);
        if(bitness == 3072)
            launch<3072>(threads, iterations, argv[1]);
        if(bitness == 4096)
            launch<4096>(threads, iterations, argv[1]);
        if(bitness == 6144)
            launch<6144>(threads, iterations, argv[1]);
        if(bitness == 8192)
            launch<8192>(threads, iterations, argv[1]);
    } catch (const std::exception& e) {
        return std::printf("Execution failed: %s.\n", e.what());
    }
    return 0;
}
