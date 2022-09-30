#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include <vector>
#include <initializer_list>

namespace HydroForest {

template<typename NumericalType>
class Vector {

protected:

    std::size_t size_;
    std::vector<NumericalType> data_;

public:

    // ---======---
    // Constructors
    // ---======---

    Vector() : size_(0), data_(0) {}
    
    Vector(std::size_t size) : size_(size), data_(size) {}

    Vector(std::size_t size, NumericalType* dataArray) : size_(size) {
        data_.assign(dataArray, dataArray + size);
    }

    Vector(std::size_t size, NumericalType value) : size_(size), data_(size, value) {}

    Vector(std::initializer_list<NumericalType> iList) : size_(iList.size()), data_(iList) {}

    Vector(const Vector& v) {
        std::size_t sizeCopy = v.size();
        std::vector<NumericalType> dataCopy = v.data();
        size_ = sizeCopy;
        data_ = dataCopy;
    }

    // Vector(Vector&& v) {
    //     *this = std::move(v);
    // }

    // ---=========================---
    // "Getter" and "Setter" functions
    // ---=========================---

    NumericalType getEntry(std::size_t index) {
        if (index > size_ || index < 0) {
            std::string errorMessage = "[HydroForest::Vector::getEntry] `index` is out of range:\n";
            errorMessage += "\tindex = " + std::to_string(index) + "\n";
            errorMessage += "\tsize = " + std::to_string(size_) + "\n";
            std::cerr << errorMessage << std::endl;
            throw std::out_of_range(errorMessage);
        }
        return data_[index];
    }

    NumericalType& operator[](std::size_t index) {
        if (index > size_ || index < 0) {
            std::string errorMessage = "[HydroForest::Vector::operator[]] `index` is out of range:\n";
            errorMessage += "\tindex = " + std::to_string(index) + "\n";
            errorMessage += "\tsize = " + std::to_string(size_) + "\n";
            std::cerr << errorMessage << std::endl;
            throw std::out_of_range(errorMessage);
        }
        return data_[index];
    }

    NumericalType& operator()(std::size_t index) {
        if (index > size_ || index < 0) {
            std::string errorMessage = "[HydroForest::Vector::operator()] `index` is out of range:\n";
            errorMessage += "\tindex = " + std::to_string(index) + "\n";
            errorMessage += "\tsize = " + std::to_string(size_) + "\n";
            std::cerr << errorMessage << std::endl;
            throw std::out_of_range(errorMessage);
        }
        return data_[index];
    }
    
    // ---============---
    // "Getter" functions
    // ---============---

    std::size_t size() const { return size_; }
    std::vector<NumericalType> data() const { return data_; }
    NumericalType* dataPointer() { return data_.data(); }

    Vector<NumericalType> getRange(std::size_t a, std::size_t b) {
        if (a > size_ || b > size_ || a < 0 || b < 0) {
            std::string errorMessage = "[HydroForest::Vector::getRange] `a` or `b` is outside of range of vector:\n";
            errorMessage += "\ta = " + std::to_string(a) + "\n";
            errorMessage += "\tb = " + std::to_string(b) + "\n";
            errorMessage += "\tsize = " + std::to_string(size_) + "\n";
            std::cerr << errorMessage << std::endl;
            throw std::out_of_range(errorMessage);
        }

        Vector<NumericalType> v((b - a) + 1);
        for (auto i = 0; i < v.size(); i++) {
            v(i) = a + i;
        }
        return v;
    }

    Vector<NumericalType> operator()(std::size_t a, std::size_t b) {
        if (a > size_ || b > size_ || a < 0 || b < 0) {
            std::string errorMessage = "[HydroForest::Vector::operator()] `a` or `b` is outside of range of vector:\n";
            errorMessage += "\ta = " + std::to_string(a) + "\n";
            errorMessage += "\tb = " + std::to_string(b) + "\n";
            errorMessage += "\tsize = " + std::to_string(size_) + "\n";
            std::cerr << errorMessage << std::endl;
            throw std::out_of_range(errorMessage);
        }

        Vector<NumericalType> v((b - a) + 1);
        for (auto i = 0; i < v.size(); i++) {
            v(i) = a + i;
        }
        return v;
    }

    Vector<NumericalType> getFromIndexSet(Vector<int> I) {
        if (I.size() > size_) {
            std::string errorMessage = "[HydroForest::Vector::operator()] `Size of index set `I` is greater than size of vector:\n";
            errorMessage += "\tI.size() = " + std::to_string(I.size()) + "\n";
            errorMessage += "\tsize = " + std::to_string(size_) + "\n";
            std::cerr << errorMessage << std::endl;
            throw std::out_of_range(errorMessage);
        }

        Vector<NumericalType> res(I.size());
        for (auto i = 0; i < I.size(); i++) {
            if (I[i] > size_ || I[i] < 0) {
                std::string errorMessage = "[HydroForest::Vector::operator()] Index in `I` is out of range:\n";
                errorMessage += "\ti = " + std::to_string(i) + "\n";
                errorMessage += "\tI[i] = " + std::to_string(I[i]) + "\n";
                errorMessage += "\tsize = " + std::to_string(size_) + "\n";
                std::cerr << errorMessage << std::endl;
                throw std::out_of_range(errorMessage);
            }
            res(i) = operator()(I(i));
        }
        return res;
    }

    Vector<NumericalType> operator()(Vector<int> I) {
        if (I.size() > size_) {
            std::string errorMessage = "[HydroForest::Vector::operator()] `Size of index set `I` is greater than size of vector:\n";
            errorMessage += "\tI.size() = " + std::to_string(I.size()) + "\n";
            errorMessage += "\tsize = " + std::to_string(size_) + "\n";
            std::cerr << errorMessage << std::endl;
            throw std::out_of_range(errorMessage);
        }

        Vector<NumericalType> res(I.size());
        for (auto i = 0; i < I.size(); i++) {
            if (I[i] > size_ || I[i] < 0) {
                std::string errorMessage = "[HydroForest::Vector::operator()] Index in `I` is out of range:\n";
                errorMessage += "\ti = " + std::to_string(i) + "\n";
                errorMessage += "\tI[i] = " + std::to_string(I[i]) + "\n";
                errorMessage += "\tsize = " + std::to_string(size_) + "\n";
                std::cerr << errorMessage << std::endl;
                throw std::out_of_range(errorMessage);
            }
            res(i) = operator()(I(i));
        }
        return res;
    }

    // ---=========---
    // Math operations
    // ---=========---

    Vector<NumericalType>& operator+=(const Vector<NumericalType>& rhs) {
        if (rhs.size() != size_) {
            std::string errorMessage = "[HydroForest::Vector::operator+=] Size of `rhs` is not the same of `this`:\n";
            errorMessage += "\trhs.size() = " + std::to_string(rhs.size()) + "\n";
            errorMessage += "\tsize = " + std::to_string(size_) + "\n";
            std::cerr << errorMessage << std::endl;
            throw std::invalid_argument(errorMessage);
        }

        for (auto i = 0; i < size_; i++) {
            data_[i] += rhs[i];
        }
        return *this;
    }

    Vector<NumericalType>& operator+=(const NumericalType rhs) {
        for (auto i = 0; i < size_; i++) {
            data_[i] += rhs;
        }
        return *this;
    }

    Vector<NumericalType> operator+(const Vector<NumericalType>& rhs) {
        if (rhs.size() != size_) {
            std::string errorMessage = "[HydroForest::Vector::operator+] Size of `rhs` is not the same of `this`:\n";
            errorMessage += "\trhs.size() = " + std::to_string(rhs.size()) + "\n";
            errorMessage += "\tsize = " + std::to_string(size_) + "\n";
            std::cerr << errorMessage << std::endl;
            throw std::invalid_argument(errorMessage);
        }

        return Vector<NumericalType>(*this) += rhs;
        
    }

    Vector<NumericalType> operator+(const NumericalType rhs) {
        return Vector<NumericalType>(*this) += rhs;
    }

    Vector<NumericalType>& operator-=(const Vector<NumericalType>& rhs) {
        if (rhs.size() != size_) {
            std::string errorMessage = "[HydroForest::Vector::operator-=] Size of `rhs` is not the same of `this`:\n";
            errorMessage += "\trhs.size() = " + std::to_string(rhs.size()) + "\n";
            errorMessage += "\tsize = " + std::to_string(size_) + "\n";
            std::cerr << errorMessage << std::endl;
            throw std::invalid_argument(errorMessage);
        }

        for (auto i = 0; i < size_; i++) {
            data_[i] -= rhs[i];
        }
        return *this;
    }

    Vector<NumericalType>& operator-=(const NumericalType& rhs) {
        for (auto i = 0; i < size_; i++) {
            data_[i] -= rhs;
        }
        return *this;
    }

    Vector<NumericalType> operator-(const Vector<NumericalType>& rhs) {
        if (rhs.size() != size_) {
            std::string errorMessage = "[HydroForest::Vector::operator-] Size of `rhs` is not the same of `this`:\n";
            errorMessage += "\trhs.size() = " + std::to_string(rhs.size()) + "\n";
            errorMessage += "\tsize = " + std::to_string(size_) + "\n";
            std::cerr << errorMessage << std::endl;
            throw std::invalid_argument(errorMessage);
        }

        return Vector<NumericalType>(*this) -= rhs;
        
    }

    Vector<NumericalType> operator-(const NumericalType rhs) {
        return Vector<NumericalType>(*this) -= rhs;
    }

    // Vector<NumericalType>& operator*=(const Vector<NumericalType>& rhs) {
    //     if (rhs.size() != size_) {
    //         std::string errorMessage = "[HydroForest::Vector::operator*=] Size of `rhs` is not the same of `this`:\n";
    //         errorMessage += "\trhs.size() = " + std::to_string(rhs.size()) + "\n";
    //         errorMessage += "\tsize = " + std::to_string(size_) + "\n";
    //         std::cerr << errorMessage << std::endl;
    //         throw std::invalid_argument(errorMessage);
    //     }

    //     for (auto i = 0; i < size_; i++) {
    //         data_[i] *= rhs[i];
    //     }
    //     return *this;
    // }

    Vector<NumericalType>& operator*=(const NumericalType& rhs) {
        for (auto i = 0; i < size_; i++) {
            data_[i] *= rhs;
        }
        return *this;
    }

    NumericalType operator*(Vector<NumericalType>& rhs) {
        if (rhs.size() != size_) {
            std::string errorMessage = "[HydroForest::Vector::operator*] Size of `rhs` is not the same of `this`:\n";
            errorMessage += "\trhs.size() = " + std::to_string(rhs.size()) + "\n";
            errorMessage += "\tsize = " + std::to_string(size_) + "\n";
            std::cerr << errorMessage << std::endl;
            throw std::invalid_argument(errorMessage);
        }

        NumericalType res = 0;
        for (auto i = 0; i < size_; i++) {
            res += rhs[i] * data_[i];
        }
        return res;
    }

    // Vector<NumericalType>& operator/=(const Vector<NumericalType>& rhs) {
    //     if (rhs.size() != size_) {
    //         std::string errorMessage = "[HydroForest::Vector::operator/=] Size of `rhs` is not the same of `this`:\n";
    //         errorMessage += "\trhs.size() = " + std::to_string(rhs.size()) + "\n";
    //         errorMessage += "\tsize = " + std::to_string(size_) + "\n";
    //         std::cerr << errorMessage << std::endl;
    //         throw std::invalid_argument(errorMessage);
    //     }

    //     for (auto i = 0; i < size_; i++) {
    //         data_[i] /= rhs[i];
    //     }
    //     return *this;
    // }

    Vector<NumericalType>& operator/=(const NumericalType& rhs) {
        for (auto i = 0; i < size_; i++) {
            data_[i] /= rhs;
        }
        return *this;
    }

};

Vector<int> vectorRange(std::size_t start, std::size_t end) {
    std::size_t N = (end - start) + 1;
    Vector<int> res(N);
    for (int i = 0; i < N; i++) {
        res[i] = start + i;
    }
    return res;
}
 
} // NAMESPACE : HydroForest

#endif // VECTOR_HPP_