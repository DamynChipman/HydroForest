#ifndef POLYNOMIAL_BASIS_HPP_
#define POLYNOMIAL_BASIS_HPP_

#include <vector>
#include <algorithm>
#include <math.h>

namespace HydroForest {

// class PolynomialBasisBase {

// public:

//     virtual int getOrder() = 0;
//     // virtual std::vector<double> getCoefficients() = 0;
//     virtual double operator()(const double x) = 0;
//     virtual std::vector<double> operator()(const std::vector<double> xs) = 0;
//     PolynomialBasisBase differentiate() {
//         std::vector<double> diffCoefs(order_ - 1);
//         for (int n = 0; n < order_ - 1; n++) {
//             diffCoefs[n] = coefs_[n+1] * n;
//         }
//         return {}
//     }

// private:

//     int order_;
//     std::vector<double> coefs_;

// };

class GeneralPolynomial {

public:

    GeneralPolynomial(std::vector<double> coefficients) : order(coefficients.size()), coefs(coefficients) { 

    }

    int getOrder() { return order; }
    std::vector<double> getCoefficients() { return coefs; }

    double operator()(const double x) {
        double f = 0;
        for (int n = 0; n < order; n++) {
            f += coefs[n] * pow(x, n);
        }
        return f;
    }

    std::vector<double> operator()(const std::vector<double> xs) {
        std::vector<double> fs(xs.size());
        std::transform(xs.begin(), xs.end(), fs.begin(), 
        [this](double x){
            return operator()(x);
        });
        return fs;
    }

    GeneralPolynomial differentiate() {
        if (order == 0) {
            return {{0}};
        }
        std::vector<double> diffCoefs(order - 1);
        for (int n = 0; n < order - 1; n++) {
            diffCoefs[n] = coefs[n+1] * n;
        }
        return {diffCoefs};
    }

// protected:

    int order;
    std::vector<double> coefs;

};

class LegendrePolynomial : public GeneralPolynomial {

public:

    LegendrePolynomial(int order) : GeneralPolynomial{std::vector<double>(order+1)} {
        generatingFunction();
    }


private:

    void generatingFunction() {
        // for (auto n = 0; n < order; n++) {
            if (this->order-1 == 0) {
                this->coefs[0] = 1;
            }
            else if (this->order-1 == 1) {
                this->coefs[1] = 1;
                this->coefs[0] = 0;
            }
            else if (this->order-1 == 2) {
                this->coefs[2] = 3.0/2.0;
                this->coefs[1] = 0;
                this->coefs[0] = 1.0/2.0;
            }
            else if (this->order-1 == 3) {
                this->coefs[3] = 5.0/2.0;
                this->coefs[2] = 0;
                this->coefs[1] = 3.0/2.0;
                this->coefs[0] = 0;
            }
            else {
                throw std::invalid_argument("NOT IMPLEMENTED!");
            }
        // }
    }

};



} // NAMESPACE: HydroForest

#endif // POLYNOMIAL_BASIS_HPP_