#ifndef POLYNOMIAL_HPP_
#define POLYNOMIAL_HPP_

#include <vector>
#include <utility>
#include <petscmat.h>

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
    std::vector<double> quadraturePoints;
    LagrangePolynomial(std::vector<double> quadraturePoints) : order(quadraturePoints.size()), quadraturePoints(quadraturePoints) {}

    Mat evalSamplePoints(std::vector<double> xSample) {
        Mat L_il;
        MatCreateDense(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, quadraturePoints.size(), xSample.size(), NULL, &L_il);

        for (int l = 0; l < xSample.size(); l++) {
            double x_l = xSample[l];
            for (int i = 0; i < quadraturePoints.size(); i++) {
                double x_i = quadraturePoints[i];
                double entry_il = 1.0;
                for (int j = 0; j < quadraturePoints.size(); j++) {
                    double x_j = quadraturePoints[j];
                    if (j != i) {
                        entry_il *= (x_l - x_j) / (x_i - x_j);
                    }
                }
                MatSetValue(L_il, i, l, entry_il, INSERT_VALUES);
            }
        }

        MatAssemblyBegin(L_il, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(L_il, MAT_FINAL_ASSEMBLY);
        return L_il;
    }

    std::pair<double, double> operator()(double x) {
        double L_il = 1.0;
        double dL_il = 0.0;
        for (auto& x_i : quadraturePoints) {
            L_il = 1.0;
            dL_il = 0.0;
            for (auto& x_j : quadraturePoints) {
                double prod = 1.0;
                if (x_i != x_j) {
                    for (auto& x_k : quadraturePoints) {
                        if (x_k != x_i && x_k != x_i) {
                            prod *= (x - x_k) / (x_i - x_k);
                        }
                    }
                    L_il *= (x - x_j) / (x_i - x_j);
                    dL_il += prod / (x_i - x_j);
                }
            }
        }
        return {L_il, dL_il};
    }

    std::pair<std::vector<double>, std::vector<double>> operator()(std::vector<double> x) {
        std::vector<double> L_il, dL_il;
        L_il.reserve(x.size());
        dL_il.reserve(x.size());
        for (auto& x_l : x) {
            std::pair<double, double> L = operator()(x_l);
            L_il.push_back(L.first);
            dL_il.push_back(L.second);
        }
        return {L_il, dL_il};
    }

};

};

#endif // POLYNOMIAL_HPP_