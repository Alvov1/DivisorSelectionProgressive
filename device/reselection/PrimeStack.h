#ifndef DIVISORSELECTIONPROGRESSIVE_PRIMESTACK_H
#define DIVISORSELECTIONPROGRESSIVE_PRIMESTACK_H

#ifdef __CUDACC__
#define gpu __host__ __device__
    #include <cuda/std/utility>
    #include <cuda/std/array>
#else
#define gpu
#include <utility>
#include <array>
#endif
#include <array>

template <typename Elem, std::size_t size>
class PrimeStack final {
#ifdef __CUDACC__
    template <typename T1, typename T2>
    using pair = cuda::std::pair<T1, T2>;
    using stack = cuda::std::array<pair<Elem, std::size_t>, size>;
#else
    template<typename T1, typename T2>
    using pair = std::pair<T1, T2>;
    using stack = std::array<pair<Elem, std::size_t>, size>;
#endif

    stack values {};
    std::size_t stackPosition = 0;

public:
    gpu constexpr PrimeStack() noexcept = default;

    gpu constexpr void push(Elem value, std::size_t index) noexcept {
        if(stackPosition >= size) return;
        values[stackPosition].first = value;
        values[stackPosition++].second = index;
    }

    gpu constexpr std::size_t pop() noexcept {
        if(stackPosition == 0) return 0;
        return values[--stackPosition].second;
    }

    gpu constexpr const auto& top() const noexcept {
        if(stackPosition > 0)
            return values[stackPosition - 1];
        else return values[0];
    }

    [[nodiscard]]
    gpu std::size_t getSize() const noexcept { return stackPosition; }

    [[nodiscard]]
    gpu constexpr bool isEmpty() const noexcept { return stackPosition == 0; }

    template <typename ReturnType> [[nodiscard]]
    gpu constexpr ReturnType unite() const noexcept {
        ReturnType value = 1u;
        for(std::size_t i = 0; i < stackPosition; ++i)
            value *= values[i].first;
        return value;
    }

    constexpr void introspect() const noexcept {
        for(std::size_t i = 0; i < stackPosition; ++i)
            std::cout << values[i].first << ' ';
        std::cout << '\n';
    }
};

#endif //DIVISORSELECTIONPROGRESSIVE_PRIMESTACK_H
