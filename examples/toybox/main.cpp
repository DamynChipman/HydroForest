#include <iostream>
#include <petsc.h>
#include <p4est.h>
#include <matplotlibcpp.h>
#include <HydroForestApp.hpp>
#include <Grid1D.hpp>
#include <Polynomial.hpp>
#include <NewtonRaphsonSolver.hpp>

namespace plt = matplotlibcpp;

struct SolveTester : public HydroForest::NewtonRaphsonSolver<SolveTester> {

    SolveTester() : NewtonRaphsonSolver<SolveTester>(*this) {}

    double objectiveFunction(double x) {
        return (4 - x) * (9 - x);
    }

    double objectiveDerivative(double x) {
        return (2*x - 13);
    }

};

int main(int argc, char** argv) {
    std::cout << "-================-" << std::endl;
    std::cout << "--=== TOYBOX ===--" << std::endl;
    std::cout << "-================-" << std::endl;

    HydroForest::HydroForestApp app(&argc, &argv);

    SolveTester solver;
    double root1 = solver.solve(3.5);
    double root2 = solver.solve(9.5);
    std::cout << "root1 = " << root1 << std::endl;
    std::cout << "root2 = " << root2 << std::endl;

    int nPlot = 400;
    HydroForest::UniformGrid1D<double> plotGrid(-1.0, 1.0, nPlot);

    // HydroForest::LegendrePolynomial legendrePoly(3);
    // std::vector<double> pLeg = legendrePoly(plotGrid.getPoints());

    // plt::plot(plotGrid.getPoints(), pLeg);
    // plt::grid(true);
    // plt::show();

    int basisOrder = 16;
    HydroForest::LobattoGrid1D<double> basisPoints(basisOrder);
    HydroForest::LagrangePolynomial L(basisPoints.getPoints());

    Mat L_il = L.evalSamplePoints(plotGrid.getPoints());

    std::vector<std::vector<double>> lagrangePolys;
    for (int i = basisOrder; i >= 0; i--) {
    // for (int i = 0; i < L.order; i++) {
        std::vector<double> L_i(nPlot);

        std::vector<int> iRows = {i};
        std::vector<int> jCols(nPlot);
        std::iota(jCols.begin(), jCols.end(), 0); // 0:nPlot
        MatGetValues(L_il, iRows.size(), iRows.data(), jCols.size(), jCols.data(), L_i.data());

        // lagrangePolys[i] = L_i;
        lagrangePolys.push_back(L_i);
    }

    for (int i = 0; i < basisOrder+1; i++) {
        std::vector<double> L_i = lagrangePolys[i];
        plt::plot(plotGrid.getPoints(), L_i, {{"label", std::to_string(i)}});
    }
    plt::plot(basisPoints.getPoints(), std::vector<double>(basisPoints.getNPoints(), 0), ".k");
    plt::title("Lagrange Polynomial Basis Functions - P = " + std::to_string(basisOrder));
    // plt::legend();
    plt::show();

    return EXIT_SUCCESS;
}