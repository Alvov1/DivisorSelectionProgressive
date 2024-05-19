#include <iostream>
#include <fstream>
#include <filesystem>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "I-AM-REALLY-STRANGE-AEU.h"
#include "Timer.h"
#include "PrimeStack.h"

using primeType = unsigned;
thrust::host_vector<primeType> loadPrimes(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load prime table: bad input file");

    const auto primesCount = (std::filesystem::file_size(fromLocation) / sizeof(primeType)) / 3;
    std::ifstream input(fromLocation, std::ios::binary);

    thrust::host_vector<primeType> primes (primesCount);
    for(auto& prime: primes)
        input.read(reinterpret_cast<char*>(&prime), sizeof(primeType));

    return primes;
}








using ValuePrecision = Aeu<512>;

__global__
void kernel(ValuePrecision* const numberAndFactor,
            const primeType* const primes,
            const std::size_t primesCount) {
    const std::size_t threadId = blockDim.x * blockIdx.x + threadIdx.x,
            threadsCount = gridDim.x * blockDim.x;






    constexpr ValuePrecision a = ValuePrecision(2u) * 1u;
    const ValuePrecision n = numberAndFactor[0];
    ValuePrecision& factor = numberAndFactor[1];




    /* Using base with some smallest primes: 2^9 * 3^6 * 5^5 * 7^4 * 11^2 * 13^2 *
      * 17 * 19 * 21 * 23 * 27 * 29 * 31 * 37. */
    constexpr ValuePrecision startingBase =
            ValuePrecision(512u) * 729u
            * 3125u * 2401u * 121u * 169u
            * 17u * 19u * 21u * 23u
            * 27u * 29u * 31u * 37u;
    constexpr auto stackSize = 10,
            primesPositionShift = 1;

    for(std::size_t i = threadId; i < primesCount; i += threadsCount) {
        /* Storing in stack: prime number with its position in array. */
        PrimeStack<primeType, stackSize + 2> stack {};
        std::size_t primesPosition = i;

        do {
            if (stack.getSize() >= stackSize) {
                const ValuePrecision test = stack.unite<ValuePrecision>() * startingBase * primes[i];

                const auto powm = ValuePrecision::powm(a, test, n);
                const auto decremented = powm - 1u;
                const auto gcd = ValuePrecision::gcd(decremented, n);

                if(gcd > 1u) {
                    char buffer[512] {};
                    gcd.getString<10>(buffer, 512, true, false);
                    printf("Thread %u. Found factor for power %s, a = 2.\n", threadId, buffer);

                    factor.tryAtomicSet(gcd);
                    return;
                }
            }

            if (stack.getSize() >= stackSize || primesPosition >= primesCount) {
                std::size_t checkedLast = primesCount;
                do {
                    primesPosition = stack.pop();
                    checkedLast -= primesPositionShift;
                } while (primesPosition == checkedLast);

                primesPosition += primesPositionShift;
                stack.push(primes[primesPosition], primesPosition);
                primesPosition += primesPositionShift;
            } else {
                stack.push(primes[primesPosition], primesPosition);
                primesPosition += primesPositionShift;
            }

        } while (!stack.isEmpty());
    }
}











int main(int argc, const char* const* const argv) {
    try {
        if(argc < 4)
            return std::printf("Usage: %s <number> <smooth location> <threads>", argv[0]);

        const ValuePrecision number = std::string_view(argv[1]);
        thrust::device_vector<ValuePrecision> numberAndFactor = thrust::host_vector{ number, ValuePrecision { 0u } };
        Timer::init() << "Factorizing number 0x" << std::hex << number << '.' << std::dec << Timer::endl;




        const thrust::device_vector<primeType> primes = loadPrimes(argv[2]);
        Timer::out << "Loaded table of primes with " << primes.size() << " elements." << Timer::endl;





        const auto threads = std::stoul(argv[3]), iterations = std::stoul(argv[4]);
        const auto timePoint = std::chrono::system_clock::now();
        const auto timeT = std::chrono::system_clock::to_time_t(timePoint);
        Timer::out << std::ctime(&timeT) << " Starting kernel <<<" << threads << ", " << threads << ">>>. Using bitness "
                   << ValuePrecision::getBitness() << '.' << Timer::endl;





        kernel<<<threads, threads>>>(
                thrust::raw_pointer_cast(numberAndFactor.data()),
                thrust::raw_pointer_cast(primes.data()),
                primes.size());






        const auto code = cudaDeviceSynchronize();
        if(code != cudaSuccess)
            return std::printf("Kernel launch failed: %s.\n", cudaGetErrorString(code));





        const ValuePrecision factor = numberAndFactor[1];
        if(factor != 0u && number % factor == 0u)
            Timer::out << "Kernel completed. Founded factor: 0x" << std::hex << number << '.' << Timer::endl;
        else Timer::out << "Kernel completed, but factor 0x" << std::hex << factor << " is incorrect." << Timer::endl;

    } catch (const std::exception& e) {
        return std::printf("Failed: %s.\n", e.what());
    }
}