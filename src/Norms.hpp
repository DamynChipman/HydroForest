#ifndef NORMS_HPP_
#define NORMS_HPP_

#include <string>
#include <cmath>
#include <algorithm>
#include "Vector.hpp"
#include "Matrix.hpp"

namespace HydroForest {

template<typename NumericalType>
NumericalType computeL1Norm(Vector<NumericalType>& a, Vector<NumericalType>& b) {
    if (a.size() != b.size()) {
        std::string errorMessage = "[HydroForest::computeL1Norm] Vectors `a` and `b` are not the same size:\n";
        errorMessage += "\ta.size() = " + std::to_string(a.size()) + "\n";
        errorMessage += "\tb.size() = " + std::to_string(b.size()) + "\n";
        std::cerr << errorMessage << std::endl;
        throw std::invalid_argument(errorMessage);
    }
    NumericalType res = 0;
    NumericalType numerator = 0;
    NumericalType denominator = 0;
    std::size_t N = a.size();
    for (auto i = 0; i < N; i++) {
        numerator += fabs(a[i] - b[i]);
        denominator += fabs(b[i]);
    }
    return numerator / denominator;
}

template<typename NumericalType>
NumericalType computeL2Norm(Vector<NumericalType>& a, Vector<NumericalType>& b) {
    if (a.size() != b.size()) {
        std::string errorMessage = "[HydroForest::computeL1Norm] Vectors `a` and `b` are not the same size:\n";
        errorMessage += "\ta.size() = " + std::to_string(a.size()) + "\n";
        errorMessage += "\tb.size() = " + std::to_string(b.size()) + "\n";
        std::cerr << errorMessage << std::endl;
        throw std::invalid_argument(errorMessage);
    }
    NumericalType res = 0;
    NumericalType numerator = 0;
    NumericalType denominator = 0;
    std::size_t N = a.size();
    for (auto i = 0; i < N; i++) {
        numerator += pow(a[i] - b[i], 2);
        denominator += pow(b[i], 2);
    }
    return pow(numerator / denominator, 0.5);
}

template<typename NumericalType>
NumericalType computeLInfNorm(Vector<NumericalType>& a, Vector<NumericalType>& b) {
    if (a.size() != b.size()) {
        std::string errorMessage = "[HydroForest::computeL1Norm] Vectors `a` and `b` are not the same size:\n";
        errorMessage += "\ta.size() = " + std::to_string(a.size()) + "\n";
        errorMessage += "\tb.size() = " + std::to_string(b.size()) + "\n";
        std::cerr << errorMessage << std::endl;
        throw std::invalid_argument(errorMessage);
    }
    NumericalType res = 0;
    NumericalType numerator = 0;
    NumericalType denominator = 0;
    std::size_t N = a.size();
    for (auto i = 0; i < N; i++) {
        numerator = std::max(numerator, fabs(a[i] - b[i]));
        denominator = std::max(denominator, fabs(b[i]));
    }
    return numerator / denominator;
}

} // NAMESPACE : HydroForest

#endif // NORMS_HPP_