#ifndef MEMORY_HPP_
#define MEMORY_HPP_

namespace HydroForest {

template<typename T>
struct MemoryBlock {
    T* ptr;
    std::size_t size;
    MemoryBlock() : ptr(nullptr), size(0) {}
    MemoryBlock(T* addr, std::size_t size) : ptr(addr), size(size) {}
    MemoryBlock(std::size_t size) : size(size) {
        ptr = (T*) new T[size];
    }
    ~MemoryBlock() {
        if (ptr != nullptr) {
            free(ptr);
        }
    }
    T& operator[](std::size_t i) {
        if (i < size) {
            return ptr[i];
        }
        else {
            std::cerr << "[HydroForest::MemoryBlock<T>::operator[]] Index is greater than size of memory block: i = " << i << ", size = " << size << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    T* data() { return ptr; }
    
};

} // NAMESPACE : HydroForest

#endif // MEMORY_HPP_