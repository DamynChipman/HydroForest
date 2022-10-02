#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

#include <HydroForestApp.hpp>
#include <Vector.hpp>
#include <Matrix.hpp>
#include <petsc.h>
#include <p4est.h>
#include <matplotlibcpp.h>
#include <Polynomial.hpp>
#include <Norms.hpp>
#include <Grid1D.hpp>

namespace plt = matplotlibcpp;

// Function to interpolate, differentiate, and integrate
static auto f = [](double x) { return cos((M_PI / 2.0) * x); };
static auto df = [](double x) { return -(M_PI/2.0)*sin((M_PI/2.0)*x); };
static double fIntegral = 4.0 / M_PI;

void runDemonstration() {
    // -----====================-----
    // Demo 1 : Show interpolation results
    // -----====================-----

    // Create plotting grid and sample function at points
    double xLower = -1;
    double xUpper = 1;
    int nPlot = 400;
    HydroForest::UniformGrid1D<double> plotGrid(xLower, xUpper, nPlot);
    HydroForest::Vector<double> samplePoints = plotGrid.getPoints();
    HydroForest::Vector<double> fSample(nPlot);
    for (auto i = 0; i < nPlot; i++) {
        fSample[i] = f(samplePoints[i]);
    }
    plt::plot(samplePoints.data(), fSample.data(), "-r");

    // Create basis points
    int basisOrder = 64;
    HydroForest::UniformGrid1D<double> uniformBasis(xLower, xUpper, basisOrder);
    HydroForest::ChebyshevGrid1D<double> chebyshevBasis(basisOrder);
    HydroForest::LegendreGrid1D<double> legendreBasis(basisOrder);
    HydroForest::LobattoGrid1D<double> lobattoBasis(basisOrder);

    // Create Lagrange nodal points and interpolation matrix
    HydroForest::Vector<double> nodalPoints = uniformBasis.getPoints();
    HydroForest::LagrangePolynomial lagrangePolyBasis(nodalPoints);
    HydroForest::Matrix<double> L_ik = lagrangePolyBasis(samplePoints);
    HydroForest::Matrix<double> L_ki = L_ik.T(); // Transpose for inner product

    HydroForest::Vector<double> fNodal(nodalPoints.size());
    for (auto i = 0; i < fNodal.size(); i++) {
        fNodal[i] = f(nodalPoints[i]);
    }
    HydroForest::Vector<double> fInterpolated = L_ki * fNodal; // f_i = L_ki * f_k

    plt::plot(samplePoints.data(), fInterpolated.data(), "--");
    plt::title("Interpolation with Uniform Basis Points - N = " + std::to_string(basisOrder));
    plt::xlabel("x");
    plt::ylabel("f(x)");
    plt::save("interpolation_example.pdf");
    plt::show();

    // -----====================-----
    // Demo 2 : Plot Lagrange polynomials at nodal basis points
    // -----====================-----

    // Create basis points and Lagrange polynomial and matrix
    HydroForest::LobattoGrid1D<double> nodalPointsBasis(basisOrder);
    HydroForest::LagrangePolynomial L(nodalPointsBasis.getPoints());
    HydroForest::Matrix<double> L_il = L(plotGrid.getPoints());

    // Get each of the Lagrange polynomials
    std::vector<HydroForest::Vector<double>> lagrangePolys;
    for (auto i = 0; i < basisOrder+1; i++) {
        HydroForest::Vector<double> L_i = L_il.getRow(i);
        lagrangePolys.push_back(L_i);
    }

    // Plot each of the Lagrange polynomials across whole domain/element
    for (int i = 0; i < basisOrder+1; i++) {
        HydroForest::Vector<double> L_i = lagrangePolys[i];
        plt::plot(plotGrid.getPoints().data(), L_i.data(), {{"label", std::to_string(i)}});
    }
    plt::plot(nodalPointsBasis.getPoints().data(), std::vector<double>(nodalPointsBasis.getNPoints(), 0), ".k");
    plt::title("Lagrange Polynomial Basis Functions - P = " + std::to_string(basisOrder));
    plt::legend();
    plt::save("lagrange_basis_functions.pdf");
    plt::show();

    return;
}

void runProblem1() {
    // Create plotting grid and sample function at points
    double xLower = -1;
    double xUpper = 1;
    int nPlot = 50;
    HydroForest::UniformGrid1D<double> plotGrid(xLower, xUpper, nPlot);
    HydroForest::Vector<double> samplePoints = plotGrid.getPoints();
    HydroForest::Vector<double> fSample(nPlot);
    HydroForest::Vector<double> dfSample(nPlot);
    for (auto i = 0; i < nPlot; i++) {
        fSample[i] = f(samplePoints[i]);
    }

    // Iterate through basis order
    HydroForest::Vector<double> uniformErrors(64);
    HydroForest::Vector<double> chebyshevErrors(64);
    HydroForest::Vector<double> legendreErrors(64);
    HydroForest::Vector<double> lobattoErrors(64);
    std::vector<HydroForest::Vector<double>*> errorVectors = {
        &uniformErrors,
        &chebyshevErrors,
        &legendreErrors,
        &lobattoErrors
    };
    HydroForest::Vector<int> basisOrders = HydroForest::vectorRange(1, 64);
    for (auto n = 0; n < basisOrders.size(); n++) {
        auto basisOrder = basisOrders[n];

        // Create basis points
        HydroForest::UniformGrid1D<double> uniformBasis(xLower, xUpper, basisOrder);
        HydroForest::ChebyshevGrid1D<double> chebyshevBasis(basisOrder);
        HydroForest::LegendreGrid1D<double> legendreBasis(basisOrder);
        HydroForest::LobattoGrid1D<double> lobattoBasis(basisOrder);
        std::vector<HydroForest::Grid1DBase<double>*> basisPoints(4);
        basisPoints[0] = &uniformBasis;
        basisPoints[1] = &chebyshevBasis;
        basisPoints[2] = &legendreBasis;
        basisPoints[3] = &lobattoBasis;

        // Iterate through basis functions
        // [uniform, chebyshev, legendre, lobatto]
        for (auto basisIndex = 0; basisIndex < basisPoints.size(); basisIndex++) {
            auto& basisPointGrid = *basisPoints[basisIndex];
            auto& errorVector = *errorVectors[basisIndex];

            // Create Lagrange nodal points and interpolation matrix
            HydroForest::Vector<double> nodalPoints = basisPointGrid.getPoints();
            HydroForest::LagrangePolynomial lagrangePolyBasis(nodalPoints);
            HydroForest::Matrix<double> L_ik = lagrangePolyBasis(samplePoints);
            HydroForest::Matrix<double> L_ki = L_ik.T(); // Transpose for inner product

            HydroForest::Vector<double> fNodal(nodalPoints.size());
            for (auto i = 0; i < fNodal.size(); i++) {
                fNodal[i] = f(nodalPoints[i]);
            }
            HydroForest::Vector<double> fInterpolated = L_ki * fNodal; // f_i = L_ki * f_k

            double l2Error = HydroForest::computeL2Norm(fInterpolated, fSample);
            errorVector[n] = l2Error;
        }
    }

    // Plot polynomial order vs norm for all basis functions
    plt::semilogy(basisOrders.data(), uniformErrors.data(), "-r");
    plt::semilogy(basisOrders.data(), chebyshevErrors.data(), "-g");
    plt::semilogy(basisOrders.data(), legendreErrors.data(), "-b");
    plt::semilogy(basisOrders.data(), lobattoErrors.data(), "-y");
    plt::title("$L_2$ Norms for Various Nodal Points\nRed : Uniform, Green : Chebyshev, Blue : Legendre, Yellow : Lobatto");
    plt::xlabel("Polynomial Order");
    plt::ylabel("Error");
    plt::save("poly_order_vs_error_function.pdf");
    plt::show();

    return;
}

void runProblem2() {
    // Create plotting grid and sample function at points
    double xLower = -1;
    double xUpper = 1;
    int nPlot = 50;
    HydroForest::UniformGrid1D<double> plotGrid(xLower, xUpper, nPlot);
    HydroForest::Vector<double> samplePoints = plotGrid.getPoints();
    HydroForest::Vector<double> fSample(nPlot);
    HydroForest::Vector<double> dfSample(nPlot);
    for (auto i = 0; i < nPlot; i++) {
        fSample[i] = f(samplePoints[i]);
        dfSample[i] = df(samplePoints[i]);
    }

    // Iterate through basis order
    HydroForest::Vector<double> uniformErrors(64);
    HydroForest::Vector<double> chebyshevErrors(64);
    HydroForest::Vector<double> legendreErrors(64);
    HydroForest::Vector<double> lobattoErrors(64);
    std::vector<HydroForest::Vector<double>*> errorVectors = {
        &uniformErrors,
        &chebyshevErrors,
        &legendreErrors,
        &lobattoErrors
    };
    HydroForest::Vector<int> basisOrders = HydroForest::vectorRange(1, 64);
    for (auto n = 0; n < basisOrders.size(); n++) {
        auto basisOrder = basisOrders[n];

        // Create basis points
        HydroForest::UniformGrid1D<double> uniformBasis(xLower, xUpper, basisOrder);
        HydroForest::ChebyshevGrid1D<double> chebyshevBasis(basisOrder);
        HydroForest::LegendreGrid1D<double> legendreBasis(basisOrder);
        HydroForest::LobattoGrid1D<double> lobattoBasis(basisOrder);
        std::vector<HydroForest::Grid1DBase<double>*> basisPoints(4);
        basisPoints[0] = &uniformBasis;
        basisPoints[1] = &chebyshevBasis;
        basisPoints[2] = &legendreBasis;
        basisPoints[3] = &lobattoBasis;

        // Iterate through basis functions
        // [uniform, chebyshev, legendre, lobatto]
        for (auto basisIndex = 0; basisIndex < basisPoints.size(); basisIndex++) {
            auto& basisPointGrid = *basisPoints[basisIndex];
            auto& errorVector = *errorVectors[basisIndex];

            // Create Lagrange nodal points and interpolation matrix
            HydroForest::Vector<double> nodalPoints = basisPointGrid.getPoints();
            HydroForest::LagrangePolynomial lagrangePolyBasis(nodalPoints);
            HydroForest::Matrix<double> dL_ik = lagrangePolyBasis.derivative(samplePoints);
            HydroForest::Matrix<double> dL_ki = dL_ik.T(); // Transpose for inner product

            HydroForest::Vector<double> fNodal(nodalPoints.size());
            for (auto i = 0; i < fNodal.size(); i++) {
                fNodal[i] = f(nodalPoints[i]);
            }
            HydroForest::Vector<double> dfInterpolated = dL_ki * fNodal; // df_i = dL_ki * f_k

            double l2Error = HydroForest::computeL2Norm(dfInterpolated, dfSample);
            errorVector[n] = l2Error;
        }
    }

    // Plot polynomial order vs norm for all basis functions
    plt::semilogy(basisOrders.data(), uniformErrors.data(), "-r");
    plt::semilogy(basisOrders.data(), chebyshevErrors.data(), "-g");
    plt::semilogy(basisOrders.data(), legendreErrors.data(), "-b");
    plt::semilogy(basisOrders.data(), lobattoErrors.data(), "-y");
    plt::title("$L_2$ Error for Function Derivative\nRed : Uniform, Green : Chebyshev, Blue : Legendre, Yellow : Lobatto");
    plt::xlabel("Polynomial Order");
    plt::ylabel("Error");
    plt::save("poly_order_vs_error_derivative.pdf");
    plt::show();

    return;
}

void runProblem3() {
    // Iterate through basis order
    HydroForest::Vector<double> legendreErrors(64);
    HydroForest::Vector<double> lobattoErrors(64);
    std::vector<HydroForest::Vector<double>*> errorVectors = {
        &legendreErrors,
        &lobattoErrors
    };
    HydroForest::Vector<int> basisOrders = HydroForest::vectorRange(1, 64);
    for (auto n = 0; n < basisOrders.size(); n++) {
        auto basisOrder = basisOrders[n];

        // Create basis points
        HydroForest::LegendreGrid1D<double> legendreBasis(basisOrder);
        HydroForest::LobattoGrid1D<double> lobattoBasis(basisOrder);
        std::vector<HydroForest::Grid1DBase<double>*> basisPoints(2);
        basisPoints[0] = &legendreBasis;
        basisPoints[1] = &lobattoBasis;

        // Iterate through basis functions
        // [legendre, lobatto]
        for (auto basisIndex = 0; basisIndex < basisPoints.size(); basisIndex++) {
            auto& basisPointGrid = *basisPoints[basisIndex];
            auto& errorVector = *errorVectors[basisIndex];

            // Create Lagrange nodal points and interpolation matrix
            HydroForest::Vector<double> nodalPoints = basisPointGrid.getPoints();
            HydroForest::Vector<double> nodalWeights = basisPointGrid.getWeights();
            HydroForest::LagrangePolynomial lagrangePolyBasis(nodalPoints);
            HydroForest::Matrix<double> L_ik = lagrangePolyBasis(nodalPoints);
            HydroForest::Matrix<double> L_ki = L_ik.T(); // Transpose for inner product

            HydroForest::Vector<double> fNodal(nodalPoints.size());
            for (auto i = 0; i < fNodal.size(); i++) {
                fNodal[i] = f(nodalPoints[i]);
            }
            HydroForest::Vector<double> fInterpolated = L_ki * fNodal; // f_i = L_ki * f_k
            double fIntegrated = fInterpolated * nodalWeights;
            double absoluteDiff = fabs(fIntegrated - fIntegral);
            errorVector[n] = absoluteDiff;
        }
    }

    // Plot polynomial order vs norm for all basis functions
    plt::semilogy(basisOrders.data(), legendreErrors.data(), "-b");
    plt::semilogy(basisOrders.data(), lobattoErrors.data(), "-y");
    plt::title("$L_2$ Error for Function Integral\nBlue : Legendre, Yellow : Lobatto");
    plt::xlabel("Polynomial Order");
    plt::ylabel("Error");
    plt::save("poly_order_vs_error_integral.pdf");
    plt::show();

    return;
}

int main(int argc, char** argv) {
    
    // -----====================-----
    // Math 597 - Project # 1
    // -----====================-----

    // Create HydroForest app
    HydroForest::HydroForestApp app(&argc, &argv);

    // -----====================-----
    // 0. Demonstration
    // -----====================-----
    runDemonstration();

    // -----====================-----
    // 1. Intergration
    // -----====================-----
    runProblem1();

    // -----====================-----
    // 2. Differentiation
    // -----====================-----
    runProblem2();

    // -----====================-----
    // 3. Integration
    // -----====================-----
    runProblem3();

    return EXIT_SUCCESS;
}