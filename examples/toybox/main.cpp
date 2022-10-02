#include <iostream>
#include <petsc.h>
#include <p4est.h>
#include <matplotlibcpp.h>
#include <HydroForestApp.hpp>
#include <Memory.hpp>
#include <Vector.hpp>
#include <Matrix.hpp>
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

    int nPlot = 400;
    HydroForest::UniformGrid1D<double> plotGrid(-1.0, 1.0, nPlot);

    int basisOrder = 4;
    HydroForest::LobattoGrid1D<double> basisPoints(basisOrder);
    HydroForest::LagrangePolynomial L(basisPoints.getPoints());
    HydroForest::Matrix<double> L_il = L(plotGrid.getPoints());

    std::vector<HydroForest::Vector<double>> lagrangePolys;
    // for (int i = basisOrder; i >= 0; i--) {
    for (auto i = 0; i < basisOrder+1; i++) {
        HydroForest::Vector<double> L_i = L_il.getRow(i);
        lagrangePolys.push_back(L_i);
    }

    for (int i = 0; i < basisOrder+1; i++) {
        HydroForest::Vector<double> L_i = lagrangePolys[i];
        plt::plot(plotGrid.getPoints().data(), L_i.data(), {{"label", std::to_string(i)}});
    }
    plt::plot(basisPoints.getPoints().data(), std::vector<double>(basisPoints.getNPoints(), 0), ".k");
    plt::title("Lagrange Polynomial Basis Functions - P = " + std::to_string(basisOrder));
    plt::legend();
    plt::show();

    return EXIT_SUCCESS;
}