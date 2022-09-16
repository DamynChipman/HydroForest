#include <iostream>
#include <petsc.h>
#include <p4est.h>
#include <matplotlibcpp.h>
#include <HydroForestApp.hpp>
#include <Grid1D.hpp>
#include <PolynomialBasis.hpp>

namespace plt = matplotlibcpp;

int main(int argc, char** argv) {
    std::cout << "-================-" << std::endl;
    std::cout << "--=== TOYBOX ===--" << std::endl;
    std::cout << "-================-" << std::endl;

    HydroForest::HydroForestApp app(&argc, &argv);

    HydroForest::UniformGrid1D<double> grid(-1, 1, 100);
    
    std::vector<double> coefs({2, 1, 2, 1});
    HydroForest::GeneralPolynomial p(coefs);
    HydroForest::GeneralPolynomial pPrime = p.differentiate();
    HydroForest::GeneralPolynomial pPrime2 = pPrime.differentiate();
    HydroForest::GeneralPolynomial pPrime3 = pPrime2.differentiate();

    std::vector<double> f = p(grid.getPoints());
    std::vector<double> fPrime = pPrime(grid.getPoints());
    std::vector<double> fPrime2 = pPrime2(grid.getPoints());
    std::vector<double> fPrime3 = pPrime3(grid.getPoints());

    plt::plot(grid.getPoints(), f);
    plt::plot(grid.getPoints(), fPrime);
    plt::plot(grid.getPoints(), fPrime2);
    plt::plot(grid.getPoints(), fPrime3);
    plt::show();

    HydroForest::LegendrePolynomial l3(3);
    HydroForest::LegendrePolynomial l2(2);
    HydroForest::LegendrePolynomial l1(1);
    HydroForest::LegendrePolynomial l0(0);

    std::vector<double> g3 = l3(grid.getPoints());
    std::vector<double> g2 = l2(grid.getPoints());
    std::vector<double> g1 = l1(grid.getPoints());
    std::vector<double> g0 = l0(grid.getPoints());

    plt::plot(grid.getPoints(), g3);
    plt::plot(grid.getPoints(), g2);
    plt::plot(grid.getPoints(), g1);
    plt::plot(grid.getPoints(), g0);
    plt::show();

    return EXIT_SUCCESS;
}