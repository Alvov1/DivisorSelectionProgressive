#ifndef AEU_MULTIPRECISION
#define AEU_MULTIPRECISION

#include <iostream>
#include <cassert>

#ifdef __CUDACC__
#define gpu __host__ __device__
    #include <cuda/std/utility>
    #include <cuda/std/array>
#else
#define gpu
#include <utility>
#include <array>
#endif

#if defined AESI_UNSAFE
#warning Enabled nightly mode for the library. Functions and methods input arguments are not checked for validity. Be really gentle
#endif


namespace {
    using byte = uint8_t;
    using block = uint32_t;

    constexpr std::size_t bitsInByte = 8, blockBitLength = sizeof(block) * bitsInByte;
    constexpr uint64_t blockBase = 1ULL << blockBitLength;
    constexpr block blockMax = 0xff'ff'ff'ffu;

#define ComparisonEqual 0
#define ComparisonLess 1
#define ComparisonGreater 2
#define ComparisonEquivalent 3
}

template <std::size_t bitness = 512> requires (bitness % blockBitLength == 0)
class Aeu final {
    static_assert(bitness > sizeof(uint64_t), "Use built-in types for numbers 64-bit or less.");

    static constexpr std::size_t blocksNumber = bitness / blockBitLength;

#ifdef __CUDACC__
    template <typename T1, typename T2>
    using pair = cuda::std::pair<T1, T2>;
    using blockLine = cuda::std::array<block, blocksNumber>;
#else
    template<typename T1, typename T2>
    using pair = std::pair<T1, T2>;
    using blockLine = std::array<block, blocksNumber>;
#endif

    /* -------------------------- @name Class members. ----------------------- */
    blockLine blocks;
    /* ----------------------------------------------------------------------- */

    /* ------------------------ @name Helper functions. ---------------------- */
    gpu static constexpr auto addLine(blockLine& dst, const blockLine& src) noexcept -> uint64_t {
        uint64_t carryOut = 0;
        for (std::size_t i = 0; i < blocksNumber; ++i) {
            uint64_t sum = static_cast<uint64_t>(dst[i]) + static_cast<uint64_t>(src[i]) + carryOut;
            carryOut = sum / blockBase; dst[i] = static_cast<block>(sum % blockBase);
        }
        return carryOut;
    }

    gpu static constexpr auto makeComplement(blockLine& line) noexcept -> uint64_t {
        uint64_t carryOut = 1;
        for(std::size_t i = 0; i < blocksNumber; ++i) {
            const uint64_t sum = blockBase - 1ULL + carryOut - static_cast<uint64_t>(line[i]);
            carryOut = sum / blockBase; line[i] = static_cast<block>(sum % blockBase);
        }
        return carryOut;
    }
    /* ----------------------------------------------------------------------- */

public:
    /* --------------------- @name Different constructors. ------------------- */

    gpu constexpr Aeu() noexcept = default;
    gpu constexpr Aeu(const Aeu& copy) noexcept = default;
    gpu constexpr Aeu& operator=(const Aeu& other) noexcept { blocks = other.blocks; return *this; }

    template <typename Integral> requires (std::is_unsigned_v<Integral>)
    gpu constexpr Aeu(Integral value) noexcept {
        if(value != 0) {
            uint64_t tValue = static_cast<uint64_t>(value);
            for (std::size_t i = 0; i < blocksNumber; ++i) {
                blocks[i] = static_cast<block>(tValue % blockBase);
                tValue /= blockBase;
            }
        } else
            blocks = {};
    }

    template <typename Char> requires (std::is_same_v<Char, char> || std::is_same_v<Char, wchar_t>)
    gpu constexpr Aeu(const Char* ptr, std::size_t size) noexcept : Aeu {} {
        if(size == 0) return;

        constexpr const Char* characters = [] {
            if constexpr (std::is_same_v<char, Char>) {
                return "09aAfFoObBxX";
            } else {
                return L"09aAfFoObBxX";
            }
        } ();
        std::size_t position = 0;
        const auto negative = (ptr[0] == [] {
            if constexpr (std::is_same_v<Char, char>) { return '-'; } else { return L'-'; }
        } ());
        if(negative) ++position;

        const auto base = [&ptr, &size, &position, &characters] {
            if (ptr[position] == characters[0] && size > position + 1) {
                switch (ptr[position + 1u]) {
                    case characters[8]:
                    case characters[9]:
                        position += 2u; return 2u;
                    case characters[6]:
                    case characters[7]:
                        position += 2u; return 8u;
                    case characters[10]:
                    case characters[11]:
                        position += 2u; return 16u;
                    default:
                        return 10u;
                }
            } else return 10u;
        } ();
        for(; position < size; ++position) {
            const auto digit = [&characters] (Char ch) {
                if(characters[0] <= ch && ch <= characters[1])
                    return static_cast<unsigned>(ch) - static_cast<unsigned>(characters[0]);
                if(characters[2] <= ch && ch <= characters[4])
                    return static_cast<unsigned>(ch) - static_cast<unsigned>(characters[2]) + 10;
                if(characters[3] <= ch && ch <= characters[5])
                    return static_cast<unsigned>(ch) - static_cast<unsigned>(characters[3]) + 10;
                return 99u;
            } (ptr[position]);

            if(digit < base) {
                this->operator*=(base);
                this->operator+=(digit);
            }
        }

        if(negative && !isZero())
            makeComplement(blocks);
    }

    template <typename Char, std::size_t arrayLength> requires (arrayLength > 1 && (std::is_same_v<Char, char> || std::is_same_v<Char, wchar_t>))
    gpu constexpr Aeu(const Char (&literal)[arrayLength]) noexcept : Aeu(literal, arrayLength) {}

    template <typename String, typename Char = typename String::value_type> requires (std::is_same_v<std::basic_string<Char>,
            typename std::decay<String>::type> || std::is_same_v<std::basic_string_view<Char>, typename std::decay<String>::type>)
    gpu constexpr Aeu(String&& stringView) noexcept : Aeu(stringView.data(), stringView.size()) {}

    /* ----------------------------------------------------------------------- */


    /* --------------------- @name Arithmetic operators. --------------------- */
    /* ------------------------- @name Unary operators. -------------------------- */
    [[nodiscard]]
    gpu constexpr auto operator+() const noexcept -> Aeu { return *this; }

    [[nodiscard]]
    gpu constexpr auto operator-() const noexcept -> Aeu {
        Aeu result = *this; makeComplement(result.blocks); return result;
    }

    gpu constexpr auto operator++() noexcept -> Aeu& { return this->operator+=(1); }

    gpu constexpr auto operator++(int) & noexcept -> Aeu {
        Aeu old = *this; operator++(); return old;
    }

    gpu constexpr auto operator--() noexcept -> Aeu& { return this->operator-=(1); }

    gpu constexpr auto operator--(int) & noexcept -> Aeu {
        Aeu old = *this; operator--(); return old;
    }
    /* --------------------------------------------------------------------------- */


    /* ------------------------ @name Addition operators. ------------------------ */
    [[nodiscard]]
    gpu constexpr auto operator+(const Aeu& addendum) const noexcept -> Aeu {
        Aeu result = *this; result += addendum; return result;
    }

    gpu constexpr auto operator+=(const Aeu& addendum) noexcept -> Aeu& {
        addLine(blocks, addendum.blocks); return *this;
    }
    /* --------------------------------------------------------------------------- */


    /* ----------------------- @name Subtraction operators. ---------------------- */
    [[nodiscard]]
    gpu constexpr auto operator-(const Aeu& subtrahend) const noexcept -> Aeu {
        Aeu result = *this; result -= subtrahend; return result;
    }

    gpu constexpr auto operator-=(const Aeu& subtrahend) noexcept -> Aeu& {
        return this->operator+=(-subtrahend);
    }
    /* --------------------------------------------------------------------------- */


    /* --------------------- @name Multiplication operators. --------------------- */
    [[nodiscard]]
    gpu constexpr auto operator*(const Aeu& factor) const noexcept -> Aeu {
        Aeu result = *this; result *= factor; return result;
    }

    gpu constexpr auto operator*=(const Aeu& factor) noexcept -> Aeu& {
        constexpr auto multiplyLines =
                [] (const blockLine& longerLine, const std::size_t longerLength,
                                           const blockLine& smallerLine, const std::size_t smallerLength) {
            blockLine buffer {};

            for(std::size_t i = 0; i < longerLength; ++i) {
                uint64_t tBlock = longerLine[i], carryOut = 0;

                for(std::size_t j = 0; j < smallerLength && i + j < buffer.size(); ++j) {
                    const auto product = tBlock * static_cast<uint64_t>(smallerLine[j]) + carryOut;
                    const auto block = static_cast<uint64_t>(buffer[i + j]) + (product % blockBase);
                    carryOut = product / blockBase + block / blockBase;
                    buffer[i + j] = block % blockBase;
                }

                if(smallerLength < blocksNumber && smallerLength + i < buffer.size())
                    buffer[smallerLength + i] += carryOut;
            }

            return buffer;
        };

        const std::size_t thisLength = this->filledBlocksNumber(), valueLength = factor.filledBlocksNumber();
        if(thisLength > valueLength)
            blocks = multiplyLines(blocks, thisLength, factor.blocks, valueLength);
        else
            blocks = multiplyLines(factor.blocks, valueLength, blocks, thisLength);

        return *this;
    }
    /* --------------------------------------------------------------------------- */


    /* ------------------------ @name Division operators. ------------------------ */
    [[nodiscard]]
    gpu constexpr auto operator/(const Aeu& divisor) const noexcept -> Aeu {
        Aeu quotient, _; divide(*this, divisor, quotient, _); return quotient;
    }

    gpu constexpr auto operator/=(const Aeu& divisor) noexcept -> Aeu& {
        return this->operator=(this->operator/(divisor));
    }
    /* --------------------------------------------------------------------------- */


    /* ------------------------- @name Modulo operators. ------------------------- */
    [[nodiscard]]
    gpu constexpr auto operator%(const Aeu& modulo) const noexcept -> Aeu {
        Aeu _, remainder; divide(*this, modulo, _, remainder); return remainder;
    }

    gpu constexpr auto operator%=(const Aeu& modulo) noexcept -> Aeu& {
        return this->operator=(this->operator%(modulo));
    }
    /* --------------------------------------------------------------------------- */
    /* ----------------------------------------------------------------------- */


    /* ----------------------- @name Bitwise operators. ---------------------- */
    template <typename Unsigned> requires (std::is_integral_v<Unsigned> && std::is_unsigned_v<Unsigned>) [[nodiscard]]
    gpu constexpr auto operator<<(Unsigned bitShift) const noexcept -> Aeu {
        Aeu result = *this; result.operator<<=(bitShift); return result;
    }

    template <typename Unsigned> requires (std::is_integral_v<Unsigned> && std::is_unsigned_v<Unsigned>)
    gpu constexpr auto operator<<=(Unsigned bitShift) noexcept -> Aeu& {
        if(bitShift >= bitness || bitShift == 0) return *this;

        const std::size_t quotient = bitShift / blockBitLength, remainder = bitShift % blockBitLength;
        const block stamp = (1u << (blockBitLength - remainder)) - 1;

        for (long long i = blocksNumber - 1; i >= (quotient + (remainder ? 1 : 0)); --i) {
            if(i - quotient > blocksNumber || i - quotient - (remainder ? 1 : 0) > blocksNumber) continue;

            blocks[i] = ((blocks[i - quotient] & stamp) << remainder) |
                        ((blocks[i - quotient - (remainder ? 1 : 0)] & ~stamp)
                                >> ((blockBitLength - remainder) % blockBitLength));
        }

        if(quotient < blocksNumber)
            blocks[quotient] = (blocks[0] & stamp) << remainder;

        for (std::size_t i = 0; i < quotient && i < blocksNumber; ++i)
            blocks[i] = 0;
        return *this;
    }
    /* ----------------------------------------------------------------------- */


    /* --------------------- @name Comparison operators. --------------------- */
    /* ------------------------ @name Equality operators. ------------------------ */
    gpu constexpr auto operator==(const Aeu& other) const noexcept -> bool { return blocks == other.blocks; };
    /* --------------------------------------------------------------------------- */


    /* ----------------------- @name Comparison operators. ----------------------- */
    gpu constexpr auto compareTo(const Aeu& other) const noexcept -> unsigned {
        for(long long i = blocksNumber - 1; i >= 0; --i) {
            const block thisBlock = blocks[i];
            const block otherBlock = other.blocks[i];

            if(thisBlock > otherBlock)
                return ComparisonGreater;

            if(thisBlock < otherBlock)
                return ComparisonLess;
        }
        return ComparisonEqual;
    };
    /* --------------------------------------------------------------------------- */


    /* ------------------------ @name Spaceship operators. ----------------------- */
    gpu constexpr auto operator!=(const Aeu& value) const noexcept -> bool {
        return !this->operator==(value);
    }

    gpu constexpr auto operator<(const Aeu& value) const noexcept -> bool {
        return this->compareTo(value) == ComparisonLess;
    }

    gpu constexpr auto operator<=(const Aeu& value) const noexcept -> bool {
        return !this->operator>(value);
    }

    gpu constexpr auto operator>(const Aeu& value) const noexcept -> bool {
        return this->compareTo(value) == ComparisonGreater;
    }

    gpu constexpr auto operator>=(const Aeu& value) const noexcept -> bool {
        return !this->operator<(value);
    }
    /* --------------------------------------------------------------------------- */
    /* ----------------------------------------------------------------------- */


    /* ---------------------- @name Supporting methods. ---------------------- */

    gpu constexpr auto setBit(std::size_t index, bool bit) noexcept -> void {
        if(index >= bitness) return; //index < bitness (bitness = blocksNumber * blockBitLength);

        const std::size_t blockNumber = index / blockBitLength, bitNumber = index % blockBitLength;
//        assert(blockNumber < blocksNumber && bitNumber < blockBitLength);

        if(bit) {
            blocks[blockNumber] |= (1u << bitNumber);
        }
        else
            blocks[blockNumber] &= (~(1u << bitNumber));
    }


    [[nodiscard]]
    gpu constexpr auto getBit(std::size_t index) const noexcept -> bool {
        if(index >= bitness) return false;
        const std::size_t blockNumber = index / blockBitLength, bitNumber = index % blockBitLength;

//        assert(blockNumber < blocksNumber && bitNumber < blockBitLength);
        return blocks[blockNumber] & (1u << bitNumber);
    }

    [[nodiscard]]
    gpu constexpr auto isZero() const noexcept -> bool { return filledBlocksNumber() == 0; }

    [[nodiscard]]
    gpu constexpr auto filledBlocksNumber() const noexcept -> std::size_t {
        for(long long i = blocksNumber - 1; i >= 0; --i)
            if(blocks[i]) return i + 1;
        return 0;
    }

    [[nodiscard]]
    gpu static constexpr auto getBitness() noexcept -> std::size_t { return bitness; }
    /* ----------------------------------------------------------------------- */


    /* -------------- @name Public arithmetic and number theory. ------------- */
    gpu constexpr auto lshiftOnce() noexcept -> void {
        block shift_out = 0;

        for (std::size_t i = 0; i < blocksNumber; i++) {
            block d = blocks[i];

            blocks[i] = shift_out | (d << 1u);
            shift_out = 1u & (d >> (blockBitLength - 1u));
        }
    }


    gpu static constexpr auto divide(const Aeu& number, const Aeu& divisor, Aeu& quotient, Aeu& remainder) noexcept -> void {
        const auto ratio = number.compareTo(divisor);

        // quotient = Aeu {}; remainder = Aeu {};

        if(ratio == ComparisonGreater) {
            // const auto bitsUsed = number.filledBlocksNumber() * blockBitLength;
            // for(long long i = bitsUsed - 1; i >= 0; --i) {
            //     remainder.lshiftOnce();
            //     bool bit = number.getBit(i);
            //     remainder.setBit(0, bit);

            //     if(remainder >= divisor) {
            //         remainder -= divisor;
            //         quotient.setBit(i, true);
            //     }
            // }
        } else if(ratio == ComparisonLess)
            remainder = number; else quotient = 1u;
    }

    gpu static constexpr auto gcd(const Aeu& first, const Aeu& second) noexcept -> Aeu {
        Aeu gcd, tGcd, quotient, remainder;

        const auto ratio = first.compareTo(second);
        if(ratio == ComparisonGreater) {
            gcd = second;
            divide(first, second, quotient, remainder);
        } else {
            gcd = first;
            divide(second, first, quotient, remainder);
        }

        for(Aeu tX = 1u, tY = 0u; remainder != 0u; ) {
            tGcd = gcd; gcd = remainder;
            divide(tGcd, gcd, quotient, remainder);
        }

        return gcd;
    }

    [[nodiscard]]
    gpu static constexpr auto powm(const Aeu& base, const Aeu& power, const Aeu& mod) noexcept -> Aeu {
        Aeu output = 1u, _, b;
        divide(base, mod, _, b);

        for(std::size_t iteration = 0u; power.filledBlocksNumber() * blockBitLength != iteration; iteration++) {
            if(power.getBit(iteration)) {
                Aeu _2, remainder;
                divide(output * b, mod, _2, remainder);
                output = remainder;
            }

            Aeu quotient, remainder;
            divide(b * b, mod, quotient, remainder);
            b = remainder;
        }

        return output;
    }

    [[nodiscard]]
    gpu static constexpr auto power2(std::size_t e) noexcept -> Aeu { Aeu result {}; result.setBit(e, true); return result; }

    /* ----------------------------------------------------------------------- */


    /* ----------------- @name Public input-output operators. ---------------- */
    template <byte base, typename Char> requires (std::is_same_v<Char, char> || std::is_same_v<Char, wchar_t> && (base == 2 || base == 8 || base == 10 || base == 16))
    gpu constexpr auto getString(Char* const buffer, std::size_t bufferSize, bool showBase = false, bool hexUppercase = false) const noexcept -> std::size_t {
        return 0;
        if(bufferSize < 2) return 0;

        std::size_t position = 0;

        if (showBase && bufferSize > 3) {
            if constexpr (base == 2) {
                if constexpr (std::is_same_v<Char, char>) {
                    memcpy(buffer, "0b", 2 * sizeof(Char));
                } else {
                    memcpy(buffer, L"0b", 2 * sizeof(Char));
                }
                position += 2;
            } else if constexpr (base == 8) {
                if constexpr (std::is_same_v<Char, char>) {
                    memcpy(buffer, "0o", 2 * sizeof(Char));
                } else {
                    memcpy(buffer, L"0o", 2 * sizeof(Char));
                }
                position += 2;
            } else if constexpr (base == 16) {
                if constexpr (std::is_same_v<Char, char>) {
                    memcpy(buffer, "0x", 2 * sizeof(Char));
                } else {
                    memcpy(buffer, L"0x", 2 * sizeof(Char));
                }
                position += 2;
            }
        }

        if(isZero()) {
            buffer[position++] = [] { if constexpr (std::is_same_v<Char, char>) { return '0'; } else { return L'0'; }}();
            return position;
        }

        if constexpr (base == 16) {
            long long iter = blocks.size() - 1;
            for (; blocks[iter] == 0 && iter >= 0; --iter);

            if constexpr (std::is_same_v<Char, char>) {
                position += snprintf(buffer + position, bufferSize - position, (hexUppercase ? "%X" : "%x"), blocks[iter--]);
                for (; iter >= 0; --iter)
                    position += snprintf(buffer + position, bufferSize - position, (hexUppercase ? "%08X" : "%08x"), blocks[iter]);
            } else {
                position += swprintf(buffer + position, bufferSize - position, (hexUppercase ? L"%X" : L"%x"), blocks[iter--]);
                for (; iter >= 0; --iter)
                    position += swprintf(buffer + position, bufferSize - position, (hexUppercase ? L"%08X" : L"%08x"), blocks[iter]);
            }
        } else {
            const auto startPosition = position;

            Aeu copy = *this;
            while (copy != 0u && position < bufferSize) {
                Aeu quotient, remainder;
                divide(copy, base, quotient, remainder);
                if constexpr (std::is_same_v<Char, char>) {
                    buffer[position++] = '0' + remainder.template integralCast<byte>();
                } else {
                    buffer[position++] = L'0' + remainder.template integralCast<byte>();
                }
                copy = quotient;
            }

            for (std::size_t i = startPosition; i * 2 < position; ++i) {
                Char t = buffer[i]; buffer[i] = buffer[position - 1 - i]; buffer[position - 1 - i] = t;
            }
        }

        buffer[position++] = Char {};
        return position;
    }

    template <typename Char> requires (std::is_same_v<Char, char> || std::is_same_v<Char, wchar_t>)
    friend constexpr auto operator<<(std::basic_ostream<Char>& ss, const Aeu& value) -> std::basic_ostream<Char>& {
        auto flags = ss.flags();

        const auto base = [] (long baseField, std::basic_ostream<Char>& ss, bool showbase) {
            unsigned base = (baseField == std::ios::hex ? 16u : (baseField == std::ios::oct ? 8u : 10u));
            if(showbase && base != 10u)
                ss << [&base] { if constexpr (std::is_same_v<Char, char>) { return base == 8 ? "0o" : "0x"; } else { return base == 8 ? L"0o" : L"0x"; }} () << std::noshowbase ;
            return base;
        } (flags & std::ios::basefield, ss, flags & std::ios::showbase);

        if(value.isZero())
            return ss << '0';

        if(base == 16u) {
            long long iter = value.blocks.size() - 1;
            for(; value.blocks[iter] == 0 && iter >= 0; --iter);

            ss << value.blocks[iter--];
            for (; iter >= 0; --iter) {
                ss.fill([] { if constexpr (std::is_same_v<Char, char>) { return '0'; } else { return L'0'; } } ());
                ss.width(8); ss << std::right << value.blocks[iter];
            }
        } else {
//            /* Well, here we use a pre-calculated magic number to ratio the length of numbers in decimal or octal notation according to bitness.
//             * * It is 2.95-98 for octal and 3.2 for decimal. */
//            constexpr auto bufferSize = static_cast<std::size_t>(static_cast<double>(bitness) / 2.95);
//            Char buffer [bufferSize] {}; std::size_t filled = 0;
//
//            Aeu copy = value;
//            while(!copy.isZero() && filled < bufferSize) {
//                Aeu quotient, remainder;
//                divide(copy, base, quotient, remainder);
//                buffer[filled++] = [] { if constexpr (std::is_same_v<Char, char>) { return '0'; } else { return L'0'; } } () + remainder.template integralCast<byte>();
//                copy = quotient;
//            }
//
//            for(; filled > 0; --filled)
//                ss << buffer[filled - 1];
        }

        return ss;
    }

    gpu constexpr void introspect() const noexcept {
        for(std::size_t i = 0; i < blocksNumber; ++i) {
            printf("%u ", blocks[i]);
        }
        printf("\n");
    }


    template <typename Char> requires (std::is_same_v<Char, char> || std::is_same_v<Char, wchar_t>)
    constexpr auto readBinary(std::basic_istream<Char>& istream, bool bigEndian = true) -> void {
        blocks = {};
        if(bigEndian) {
            for(auto it = blocks.rbegin(); it != blocks.rend(); ++it)
                if(!istream.read(reinterpret_cast<char*>(&*it), sizeof(block))) break;
        } else {
            for(auto& tBlock: blocks)
                if(!istream.read(reinterpret_cast<char*>(&tBlock), sizeof(block))) break;
        }
    }
    /* ----------------------------------------------------------------------- */


    /* -------------------- @name Public casting operators. ------------------ */
    template <typename Integral> requires (std::is_integral_v<Integral>) [[nodiscard]]
    gpu constexpr auto integralCast() const noexcept -> Integral {
        const uint64_t value = (static_cast<uint64_t>(blocks[1]) << blockBitLength) | static_cast<uint64_t>(blocks[0]);
        return static_cast<Integral>(value % (sizeof(Integral) * 8));
    }

    template <std::size_t newBitness> requires (newBitness != bitness) [[nodiscard]]
    gpu constexpr auto precisionCast() const noexcept -> Aeu<newBitness> {
        Aeu<newBitness> result {};

        const std::size_t blockBoarder = (newBitness > bitness ?
                Aeu<bitness>::blocksNumber
                :
                Aeu<newBitness>::totalBlocksNumber());

        for(std::size_t blockIdx = 0; blockIdx < blockBoarder; ++blockIdx)
            result.setBlock(blockIdx, blocks[blockIdx]);

        return result;
    }
    /* ----------------------------------------------------------------------- */


#if defined __CUDACC__
    /* ------------------- @name Atomic-like CUDA operators. ----------------- */
    __device__ constexpr auto tryAtomicSet(const Aeu& value) noexcept -> void {
        for(std::size_t i = 0; i < blocksNumber; ++i)
            atomicExch(&blocks[i], value.blocks[i]);
    }
    /* ----------------------------------------------------------------------- */
#endif
};

#endif //AEU_MULTIPRECISION
