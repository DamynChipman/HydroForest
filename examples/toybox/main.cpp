#include <iostream>
#include <petsc.h>
#include <p4est.h>
#include <matplotlibcpp.h>
#include <HydroForestApp.hpp>
#include <Grid1D.hpp>
#include <Polynomial.hpp>

namespace plt = matplotlibcpp;

int main(int argc, char** argv) {
    std::cout << "-================-" << std::endl;
    std::cout << "--=== TOYBOX ===--" << std::endl;
    std::cout << "-================-" << std::endl;

    HydroForest::HydroForestApp app(&argc, &argv);

    int nPlot = 100;
    HydroForest::UniformGrid1D<double> plotGrid(-1, 1, nPlot);

    // HydroForest::LegendrePolynomial legendrePoly(3);
    // std::vector<double> pLeg = legendrePoly(plotGrid.getPoints());

    // plt::plot(plotGrid.getPoints(), pLeg);
    // plt::grid(true);
    // plt::show();

    int nBasisPoints = 2;
    HydroForest::LegendreGrid1D<double> basisPoints(nBasisPoints);
    HydroForest::LagrangePolynomial L(basisPoints.getPoints());

    Mat L_il = L.evalSamplePoints(plotGrid.getPoints());

    std::vector<std::vector<double>> lagrangePolys;
    // for (int i = nBasisPoints; i >= 0; i--) {
    for (int i = 0; i < L.order; i++) {
        std::vector<double> L_i(nPlot);

        std::vector<int> iRows = {i};
        std::vector<int> jCols(nPlot);
        std::iota(jCols.begin(), jCols.end(), 0); // 0:nPlot
        MatGetValues(L_il, iRows.size(), iRows.data(), jCols.size(), jCols.data(), L_i.data());

        // lagrangePolys[i] = L_i;
        lagrangePolys.push_back(L_i);
    }

    for (int i = 0; i < nBasisPoints+1; i++) {
        std::vector<double> L_i = lagrangePolys[i];
        plt::plot(plotGrid.getPoints(), L_i, {{"label", std::to_string(i)}});
    }
    plt::plot(basisPoints.getPoints(), std::vector<double>(basisPoints.getNPoints(), 0), ".k");
    plt::title("Lagrange Polynomial Basis Functions - P = " + std::to_string(nBasisPoints));
    // plt::legend();
    plt::show();

    return EXIT_SUCCESS;
}