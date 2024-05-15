#ifndef CIRCLE_BUFFER_H
#define CIRCLE_BUFFER_H

#ifdef __CUDACC__
    #define gpu __host__ __device__
    #include <cuda/std/utility>
    #include <cuda/std/array>
#else
#define gpu
    #include <utility>
    #include <array>
#endif

template <typename Elem, std::size_t length>
class CircularBuffer {
    std::array<Elem, length> buffer {};   /// pointer to allocated memory buffer
    std::size_t filledElements {};  /// number of items in the queue
    std::size_t putIndex {};    /// array index where the next push_back will occur
    std::size_t takeIndex {};   /// array index where next pop_front will occur

public:
    /**
     * Push an item into the queue
     *
     * Place the item passed into the end of the queue. The item will be returned
     * to the calling program, in FIFO order, using pop_front().
     * If the buffer is already full before the push_back(), the behavior will depend on the
     * the setting controlled by the setFullOverwrite() method.
     *
     * @param itm    a pointer to data buffer of the item to be saved. Data size must be size specified in the constructor.
     * @return true  if the item was successfully placed in the queue, false otherwise
     */
    gpu constexpr bool push_back(Elem itm) noexcept {
        if (isFull())
            Elem value = pop_front();

        // Save item and adjust the tail pointer
        buffer[putIndex++] = itm;
        filledElements++;
        if (putIndex == length)
            putIndex = 0;

        return(true);
    }

    /**
     * Pop an item from the queue
     *
     * Return the first available item in the queue, copied into the buffer specified,
     * returning a pointer to the copied item. If no items are available (queue is
     * empty), then no data is copied and the method returns a NULL pointer.
     *
     * @param itm  a pointer to data buffer for the retrieved item to be saved. Data size must be size specified in the constructor.
     * @return pointer to the memory buffer or NULL if the queue is empty
     */
    gpu Elem pop_front() {
        if (isEmpty())
            return(Elem());

        // Copy data from the buffer
        Elem value = buffer[takeIndex++];
        filledElements--;

        // If head has reached last item, wrap it back around to the start
        if (takeIndex == length)
            takeIndex = 0;

        return value;
    }

    gpu bool isEmpty() const noexcept { return filledElements == 0; };

    gpu bool isFull() const noexcept { return (filledElements != 0 && filledElements == length); };

    void introspect() const noexcept {
        for(std::size_t i = 0; i < filledElements; ++i)
            std::cout << buffer[(takeIndex + i) % length] << ' ';
        std::cout << std::endl;
    }
};

#endif  // CIRCLE_BUFFER_H