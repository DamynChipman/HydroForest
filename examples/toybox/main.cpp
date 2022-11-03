#include <cmath>
#include <iostream>
#include <string>
#include <variant>
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
#include <DGFiniteElements.hpp>
#include <Mesh1D.hpp>
#include <InitialCondition.hpp>

namespace plt = matplotlibcpp;

void advanceCG(HydroForest::ElementMesh1D<double>& mesh, double tFinal) {

    // Initialize all the things
    HydroForest::HydroForestApp& app = HydroForest::HydroForestApp::getInstance();
    HydroForest::Options& options = app.getOptions();

    // Get values from mesh
    std::size_t N = mesh.order();
    std::size_t nElements = mesh.size();
    std::size_t nDOFs = N*nElements + 1;
    std::vector<HydroForest::Element1D<double>>& elements = mesh.elements();
    HydroForest::LagrangePolynomial& poly = elements[0].polynomial();
    HydroForest::Grid1DBase<double>& grid = *elements[0].grid();
    HydroForest::Grid1DBase<double>& qGrid = *elements[0].quadratureGrid();

    // Create ID matrix
    HydroForest::CGIDMatrix ID(N, nElements);
    ID(N, nElements-1) = 0;

    // Create element matrices
    HydroForest::CGMassMatrix<double> M_ij(poly, grid, qGrid);
    HydroForest::CGDerivativeMatrix<double> D_ij(poly, grid, qGrid);

    // Create matrices
    HydroForest::CGDirectStiffnessSummationOperator<double> DSS(mesh.elements(), ID);
    HydroForest::Matrix<double> M_IJ = DSS.operateWithMetricTerm(M_ij);
    HydroForest::Matrix<double> D_IJ = DSS.operate(D_ij);
    
    // Apply BC
    std::size_t iLeft = ID(0, nElements-1);
    std::size_t iRight = ID(N, nElements-1);
    M_IJ(M_IJ.nRows()-1, M_IJ.nCols()-1) = 1;
    M_IJ(iLeft, nElements) = M_IJ(iLeft, 0);
    M_IJ(iLeft, 0) = 0;
    D_IJ(iLeft, nElements) = D_IJ(iLeft, 0);
    D_IJ(iLeft, 0) = 0;

    // Compute RHS matrix
    HydroForest::Matrix<double> Dhat = HydroForest::solve(M_IJ, D_IJ);
    Dhat *= 2.0; // Put velocity into Dhat matrix
    HydroForest::Matrix<double> nDhat = -Dhat;

    if (nElements <= 4 && N <= 2) {
        std::cout << "--=== MATRICES ===--" << std::endl;
        std::cout << "nElements = " << nElements << "  N = " << N << std::endl;
        std::cout << "M_ij = " << M_ij << std::endl;
        std::cout << "D_ij = " << D_ij << std::endl;
        std::cout << "ID = " << ID << std::endl;
        std::cout << "M_IJ = " << M_IJ << std::endl;
        std::cout << "D_IJ = " << D_IJ << std::endl;
        std::cout << "Dhat = " << Dhat << std::endl;
        std::cout << "-Dhat = " << nDhat << std::endl;
    }

    // Iterate over time
    double tStart = 0;
    int nTime = 100000;
    double time = tStart;
    HydroForest::RungeKutta3<double> timeIntegrator;
    std::vector<double> plotTimes = {.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0};
    std::size_t nPlot = 0;
    bool doPlot = false;
    for (auto n = 0; n <= nTime; n++) {
        
        // Time advance logic
        double dt = 1e-2*timeIntegrator.getMaxTimeStep(2.0, mesh[0].dx());
        double CFL = 2.0*dt/mesh[0].dx();
        if (time + dt > tFinal) {
            dt = tFinal - time;
        }
        if (time + dt >= plotTimes[nPlot]) {
            dt = plotTimes[nPlot] - time;
            doPlot = true;
        }
        time += dt;
        app.log("Timestep = %i, Time = %f, dt = %f, CFL = %f", n, time, dt, CFL);

        // Construct global solution
        HydroForest::Vector<double> q_n(nDOFs, 0);
        for (auto e = 0; e < mesh.size(); e++) {
            HydroForest::Element1D<double>& element = mesh[e];
            for (auto i = 0; i < element.size(); i++) {
                int I = ID(i,e);
                q_n[I] = element.solution()[i];
            }
        }

        HydroForest::Vector<double> q_update = 1.0*timeIntegrator.update(time, dt, q_n, -Dhat);

        // Set global solution into mesh
        for (auto e = 0; e < mesh.size(); e++) {
            HydroForest::Element1D<double>& element = mesh[e];
            for (auto i = 0; i < element.size(); i++) {
                HydroForest::Vector<int> elementMap = ID.getCol(e);
                element.solution()[i] = q_update[elementMap[i]];
            }
        }

        // Post process and plot
        if (doPlot) {
            mesh.plot("-or");
            plt::title("Time = " + std::to_string(time));
            plt::xlim(-1.0, 1.0);
            plt::ylim(-0.2, 1.2);
            plt::show();

            nPlot++;
            doPlot = false;
        }

        if (time >= tFinal) {
            break;
        }
    }

}

void advanceDG(HydroForest::ElementMesh1D<double>& mesh, double tFinal) {

    // Initialize all the things
    HydroForest::HydroForestApp& app = HydroForest::HydroForestApp::getInstance();
    HydroForest::Options& options = app.getOptions();

    // Get values from mesh
    std::size_t N = mesh.order();
    std::size_t nElements = mesh.size();
    std::size_t nDOFs = nElements*(N + 1);
    std::vector<HydroForest::Element1D<double>>& elements = mesh.elements();
    HydroForest::LagrangePolynomial& poly = elements[0].polynomial();
    HydroForest::Grid1DBase<double>& grid = *elements[0].grid();
    HydroForest::Grid1DBase<double>& qGrid = *elements[0].quadratureGrid();

    // Create ID matrix
    HydroForest::DGIDMatrix ID(N, nElements);
    // ID(N, nElements-1) = 0;

    // Create element matrices
    HydroForest::DGMassMatrix<double> M_ij(poly, grid, qGrid);
    HydroForest::DGDerivativeMatrix<double> D_ij(poly, grid, qGrid);
    HydroForest::DGFluxMatrix<double> F_ij(N);

    // Create matrices
    HydroForest::DGDirectStiffnessSummationOperator<double> DSS(mesh.elements(), ID);
    HydroForest::Matrix<double> M_IJ = DSS.operateWithMetricTerm(M_ij);
    HydroForest::Matrix<double> D_IJ = DSS.operate(D_ij);
    HydroForest::Matrix<double> F_IJ = DSS.operate(F_ij);
    HydroForest::Matrix<double> Dhat = HydroForest::solve(M_IJ, D_IJ);
    HydroForest::Matrix<double> Fhat = HydroForest::solve(M_IJ, F_IJ);

    if (nElements <= 4 && N <= 2) {
        std::cout << "--=== MATRICES ===--" << std::endl;
        std::cout << "nElements = " << nElements << "  N = " << N << std::endl;
        std::cout << "M_ij = " << M_ij << std::endl;
        std::cout << "D_ij = " << D_ij << std::endl;
        std::cout << "F_ij = " << F_ij << std::endl;
        std::cout << "ID = " << ID << std::endl;
        std::cout << "M_IJ = " << M_IJ << std::endl;
        std::cout << "D_IJ = " << D_IJ << std::endl;
        std::cout << "F_IJ = " << F_IJ << std::endl;
        std::cout << "Dhat = " << Dhat << std::endl;
        std::cout << "Fhat = " << Fhat << std::endl;
    }

    // Iterate over time
    double tStart = 0;
    int nTime = 100000;
    double time = tStart;
    HydroForest::RungeKutta3<double> timeIntegrator;
    std::vector<double> plotTimes = {.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0};
    std::size_t nPlot = 0;
    bool doPlot = false;
    for (auto n = 0; n <= nTime; n++) {
        
        // Time advance logic
        double dt = 1e-2*timeIntegrator.getMaxTimeStep(2.0, mesh[0].dx());
        double CFL = 2.0*dt/mesh[0].dx();
        if (time + dt > tFinal) {
            dt = tFinal - time;
        }
        if (time + dt >= plotTimes[nPlot]) {
            dt = plotTimes[nPlot] - time;
            doPlot = true;
        }
        time += dt;
        app.log("Timestep = %i, Time = %f, dt = %f, CFL = %f", n, time, dt, CFL);

        // Construct global solution
        HydroForest::Vector<double> q_n(nDOFs, 0);
        for (auto e = 0; e < mesh.size(); e++) {
            HydroForest::Element1D<double>& element = mesh[e];
            for (auto i = 0; i < element.size(); i++) {
                int I = ID(i,e);
                q_n[I] = element.solution()[i];
            }
        }

        // Update solution in time
        HydroForest::Vector<double> q_update = timeIntegrator.update(time, dt, q_n, [&](HydroForest::Vector<double> q_n){
            HydroForest::Vector<double> f_n = 2.0*q_n;
            HydroForest::Vector<double> f_star = HydroForest::DGGlobalRusanovFluxVector<double>(mesh.elements(), ID, q_n, f_n, 2.0);
            HydroForest::Vector<double> r_n = Dhat*f_n - Fhat*f_star;
            return r_n;
        });

        // Set global solution into mesh
        for (auto e = 0; e < mesh.size(); e++) {
            HydroForest::Element1D<double>& element = mesh[e];
            for (auto i = 0; i < element.size(); i++) {
                HydroForest::Vector<int> elementMap = ID.getCol(e);
                element.solution()[i] = q_update[elementMap[i]];
                element.flux()[i] = 2.0*q_update[elementMap[i]];
            }
        }

        // Post process and plot
        if (doPlot) {
            mesh.plot("-or");
            plt::title("Time = " + std::to_string(time));
            plt::xlim(-1.0, 1.0);
            plt::ylim(-0.2, 1.2);
            plt::show();

            nPlot++;
            doPlot = false;
        }

        if (time >= tFinal) {
            break;
        }
    }

}

int main(int argc, char** argv) {
    std::cout << "-================-" << std::endl;
    std::cout << "--=== TOYBOX ===--" << std::endl;
    std::cout << "-================-" << std::endl;

    HydroForest::HydroForestApp app(&argc, &argv);
    HydroForest::Options& options = app.getOptions();
    options["scheme"] = "DG";

    std::size_t nElements = 4;
    std::size_t N = 1;
    if (argc > 1) {
        nElements = atoi(argv[1]);
        N = atoi(argv[2]);
    }
    
    double xLower = -1;
    double xUpper = 1;
    double speed = 2.0;
    HydroForest::ElementMesh1D<double> mesh(xLower, xUpper, nElements, N);
    std::vector<HydroForest::Element1D<double>>& elements = mesh.elements();

    mesh.setInitialCondition(
        [&](double x){
            return exp(-64.0*pow(x,2));
        },
        [&](double q){
            return q*speed;
        }
    );

    // Plot initial condition
    mesh.plot("-or");
    plt::title("Initial Condition");
    plt::xlim(-1.0, 1.0);
    plt::ylim(-0.2, 1.2);
    plt::show();

    double tFinal = 1.0;

    if (std::get<std::string>(options["scheme"]) == "CG") {
        advanceCG(mesh, tFinal);
    }
    else if (std::get<std::string>(options["scheme"]) == "DG") {
        advanceDG(mesh, tFinal);
    }
    else {
        throw std::invalid_argument("Invalid Scheme; Options are `CG` or `DG`");
    }

    return EXIT_SUCCESS;
}