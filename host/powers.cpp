#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>

#include "Aesi.h"
#include "Timer.h"

const struct { unsigned x {}; unsigned y {}; } gridDim = { 64, 1 }, blockDim = { 64, 1 }, blockIdx = { 0, 0 }, threadIdx = { 0, 0 };

//3,//blockDim.x * blockIdx.x + threadIdx.x,
const unsigned threadId = 0, threadsCount = gridDim.x * blockDim.x, iterationBorder = 2;
using Uns = Aeu<2048>;

std::vector<Uns> loadPowers(const std::filesystem::path& fromLocation) {
    if(!std::filesystem::is_regular_file(fromLocation))
        throw std::invalid_argument("Failed to load powers table: bad input file");

    const auto recordsNumber = static_cast<std::size_t>(std::filesystem::file_size(fromLocation) / (static_cast<std::size_t>(Uns::getBitness() / 8)));
    std::vector<Uns> records (recordsNumber);

    std::ifstream input(fromLocation, std::ios::binary);
    for(auto& value: records)
        value.readBinary(input);

    return records;
}

void kernel(const std::vector<Uns>& powers, std::pair<Uns, Uns>& numberAndFactor) {
    auto& [n, factor] = numberAndFactor;
    const auto checkFactor = [&n, &factor] (const Uns& candidate) {
        if(candidate < 2 || candidate >= n)
            return false;
        factor = candidate; return true;
    };

    uint64_t a = threadId * iterationBorder + 2;
    for(unsigned i = threadId; i < powers.size(); i += blockDim.x) {
        for(unsigned it = 0; it < iterationBorder; ++it) {
            if(!factor.isZero())
                return;

            if(checkFactor(Uns::gcd(Uns::powm(a, powers[i], n) - 1, n)))
                return;

            a += threadsCount * iterationBorder;
        }
    }
}

int main() {
    const std::filesystem::path powers = "../../powers/powers-2048.txt";
    const Uns tNumber = "0x6508b5dc2856a2dd";

    std::pair<Uns, Uns> numberAndFactor = { tNumber, { 0 } };
    Timer::init() << "Factorizing number " << std::hex << std::showbase << numberAndFactor.first <<
                  std::dec << " (" << numberAndFactor.first.bitCount() << " bits)." << Timer::endl;

    const auto values = loadPowers(powers);
    Timer::out << "Loaded table of smooth with " << values.size() << " elements." << Timer::endl;

    kernel(values, numberAndFactor);
    if (numberAndFactor.second != 0 && numberAndFactor.first % numberAndFactor.second == 0)
        Timer::out << "Found factor: " << std::hex << std::showbase << numberAndFactor.second << Timer::endl;
    else Timer::out << "Factor is not fount or incorrect" << Timer::endl;
    return 0;
}