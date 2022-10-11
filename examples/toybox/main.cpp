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
#include <Element1D.hpp>
#include <TimeIntegration.hpp>
#include <NewtonRaphsonSolver.hpp>
#include <CGFiniteElements.hpp>

namespace plt = matplotlibcpp;

int main(int argc, char** argv) {
    std::cout << "-================-" << std::endl;
    std::cout << "--=== TOYBOX ===--" << std::endl;
    std::cout << "-================-" << std::endl;

    HydroForest::HydroForestApp app(&argc, &argv);
    
    double xLower = -1;
    double xUpper = 1;
    std::size_t nElements = 8;
    HydroForest::UniformGrid1D<double> elementGrid(xLower, xUpper, nElements);
    double dx = elementGrid[1] - elementGrid[0];

    std::size_t N = 4;
    HydroForest::LobattoGrid1D<double> grid(N);
    HydroForest::LagrangePolynomial poly(grid.getPoints());
    HydroForest::LobattoGrid1D<double> qGrid(N+1);

    HydroForest::CGMassMatrix<double> M(poly, grid, qGrid);
    // M *= 1.0 / (dx / 6.0);
    HydroForest::CGDerivativeMatrix<double> D(poly, grid, qGrid);
    HydroForest::CGIDMatrix ID(N, nElements);

    std::cout << grid << std::endl;
    std::cout << "dx = " << dx << std::endl;
    std::cout << "M = " << M << std::endl;
    std::cout << "D = " << D << std::endl;
    std::cout << "ID = " << ID << std::endl;


    std::vector<HydroForest::Element1D<double>> elements = HydroForest::createElementGrid(elementGrid, grid, qGrid);

    HydroForest::Matrix<double> M_Global = HydroForest::directStiffnessSummation<double>(elements, ID, M);
    HydroForest::Matrix<double> D_Global = HydroForest::directStiffnessSummation<double>(elements, ID, D);

    std::cout << "M_Global = " << M_Global << std::endl;
    std::cout << "D_Global = " << D_Global << std::endl;

    return EXIT_SUCCESS;
}