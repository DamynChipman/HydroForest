#ifndef POLYNOMIAL_HPP_
#define POLYNOMIAL_HPP_

#include <vector>
#include <utility>
#include <petscmat.h>

#include "Vector.hpp"
#include "Matrix.hpp"

namespace HydroForest {

/**
 * @brief Legendre polynomial construction
 * 
 * Builds an object that can be called via the `()` operator to evaluate a point `x`.
 * Legendre polynomials are formed via the recursive expression:
 *      `\phi^{Leg}_0 = 1`
 *      `\phi^{Leg}_1 = x`
 *      `\phi^{Leg}_N (x) = \frac{2N - 1}{N} x \phi^{Leg}_{N-1} (x) - \frac{N - 1}{N} \phi^{Leg}_{N-2} (x)`
 * 
 */
struct LegendrePolynomial {

    /**
     * @brief Order of the polynomial
     * 
     */
    int order = 0;

    /**
     * @brief Construct a new Legendre Polynomial object
     * 
     */
    LegendrePolynomial() {}

    /**
     * @brief Construct a new Legendre Polynomial object
     * 
     * @param order Order of the polynomial
     */
    LegendrePolynomial(int order) : order(order) {}

    /**
     * @brief Computes the Legendre polynomial, it's first and second derivatives at the given point `x`
     * 
     * @param x Point to evaluate
     * @return std::vector<double> {l0, l0_1, l0_2} where `_n` denotes the n-th derivative
     */
    double operator()(double x) {
        double l1 = 0.0;
        double l0 = 1.0;
        for (int i = 1; i <= order; i++) {
            double ii = (double) i;
            double l2 = l1;
            double a = (2.0*ii - 1.0) / ii;
            double b = (ii - 1.0) / ii;
            l1 = l0;
            l0 = a*x*l1 - b*l2;
        }
        return l0;
    }

    std::vector<double> operator()(std::vector<double> x) {
        std::vector<double> f(x.size());
        for (int l = 0; l < x.size(); l++) {
            f[l] = operator()(x[l]);
        }
        return f;
    }

    std::vector<double> evalD012(double x) {
        double l1 = 0.0; double l1_1 = 0.0; double l1_2 = 0.0;
        double l0 = 1.0; double l0_1 = 0.0; double l0_2 = 0.0;
        for (int i = 1; i <= order; i++) {
            double ii = (double) i;
            double l2 = l1; double l2_1 = l1_1; double l2_2 = l1_2;
            double a = (2.0*ii - 1.0) / ii;
            double b = (ii - 1.0) / ii;
            l1 = l0;
            l1_1 = l0_1;
            l1_2 = l0_2;
            l0 = a*x*l1 - b*l2;
            l0_1 = a*(l1 + x*l1_1) - b*(l2_1);
            l0_2 = a*(2.0*l1_1 + x*l1_2) - b*l2_2;
        }
        return {l0, l0_1, l0_2};
    }

};

struct LagrangePolynomial {

    int order;
    Vector<double> nodalPoints;
    LagrangePolynomial(Vector<double> nodalPoints) : order(nodalPoints.size()), nodalPoints(nodalPoints) {}

    Matrix<double> operator()(Vector<double> xSample) {
        std::size_t Q = xSample.size();
        std::size_t N = nodalPoints.size();
        Matrix<double> L_il(N, Q);
        for (auto l = 0; l < Q; l++) {
            double x_l = xSample[l];
            for (auto i = 0; i < N; i++) {
                double x_i = nodalPoints[i];
                L_il(i, l) = 1.0;
                for (auto j = 0; j < N; j++) {
                    if (i != j) {
                        double x_j = nodalPoints[j];
                        L_il(i,l) *= (x_l - x_j) / (x_i - x_j);
                    }
                }
            } 
        }
        return L_il;
    }

    Matrix<double> derivative(Vector<double> xSample) {
        std::size_t Q = xSample.size();
        std::size_t N = nodalPoints.size();
        Matrix<double> dL_il(N, Q);
        for (auto l = 0; l < Q; l++) {
            double x_l = xSample[l];
            for (auto i = 0; i < N; i++) {
                double x_i = nodalPoints[i];
                for (auto j = 0; j < N; j++) {
                    double x_j = nodalPoints[j];
                    double prod = 1.0;
                    if (j != i) {
                        for (auto k = 0; k < N; k++) {
                            double x_k = nodalPoints[k];
                            if (k != i && k != j) {
                                prod *= (x_l - x_k) / (x_i - x_k);
                            }
                        }
                        dL_il(i,l) += prod / (x_i - x_j);
                    }
                }
            }
        }
        return dL_il;
    }

};

};

#endif // POLYNOMIAL_HPP_