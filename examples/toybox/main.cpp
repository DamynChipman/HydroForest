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

// void run(HydroForest::Mesh1DBase<double>& mesh, HydroForest::)

// void run(HydroForest::ElementMesh1D<double>& mesh) {

//     // Initialize all the things
//     HydroForest::HydroForestApp& app = HydroForest::HydroForestApp::getInstance();
//     HydroForest::Options& options = app.getOptions();

//     // Construct M_e_ij, Dtilde_e_ij, F_e_ij
//     // HydroForest::Matrix<double> M_e_ij = HydroForest::CGMassMatrix<double>()
//     for (auto e = 0; e < mesh.size(); e++) {
//         HydroForest::Element1D<double>& element = mesh[e];
//         element.matrices()["mass-matrix"] = HydroForest::CGMassMatrix<double>(element.polynomial(), *element.grid(), *element.quadratureGrid());
//         element.matrices()["weak-derivative-matrix"] = HydroForest::DGDerivativeMatrix<double>(element.polynomial(), *element.grid(), *element.quadratureGrid());
//         element.matrices()["strong-derivative-matrix"] = HydroForest::CGDerivativeMatrix<double>(element.polynomial(), *element.grid(), *element.quadratureGrid());
//         element.matrices()["flux-matrix"] = HydroForest::DGFluxMatrix<double>(mesh.order());

//         if (e == 0) {
//             std::cout << "mass-matrix = " << element.matrices()["mass-matrix"] << std::endl;
//             std::cout << "weak-derivative-matrix = " << element.matrices()["weak-derivative-matrix"] << std::endl;
//             std::cout << "strong-derivative-matrix = " << element.matrices()["strong-derivative-matrix"] << std::endl;
//             std::cout << "flux-matrix = " << element.matrices()["flux-matrix"] << std::endl;
//         }
//     }

//     // Create ID matrices
//     HydroForest::Matrix<int> ID;
//     if (std::get<std::string>(options["scheme"]) == "CG") {
//         ID = HydroForest::CGIDMatrix(mesh.order(), mesh.nElements());
//         ID(mesh.order(), mesh.nElements()-1) = 0;
//     }
//     else if (std::get<std::string>(options["scheme"]) == "DG") {

//     }
//     else {
//         throw std::invalid_argument("blah");
//     }

//     // Construct global M_IJ
//     HydroForest::Matrix<double> M_IJ;
//     if (std::get<std::string>(options["scheme"]) == "CG") {
//         // Construct via DSS
//         M_IJ = HydroForest::directStiffnessSummationMatrix(mesh.elements(), ID, "mass-matrix");
//         M_IJ(M_IJ.nRows()-1, M_IJ.nCols()-1) = 1;
//     }
//     else if (std::get<std::string>(options["scheme"]) == "DG") {

//     }
//     else {
//         throw std::invalid_argument("blah");
//     }
//     std::cout << "M_IJ = " << M_IJ << std::endl;

//     // Initialize RHS vector
//     std::size_t nDOFs;
//     if (std::get<std::string>(options["scheme"]) == "CG") {
//         nDOFs = mesh.nElements() * mesh.order() + 1;
//     }
//     else if (std::get<std::string>(options["scheme"]) == "DG") {
//         nDOFs = mesh.nElements() * mesh.order() + mesh.nElements();
//     }
//     else {
//         throw std::invalid_argument("blah");
//     }
//     // HydroForest::Vector<double> R_I(nDOFs, 0);

//     // Iterate over time
//     double tStart = 0;
//     double tEnd = 1;
//     int nTime = 100;
//     double dt = (tEnd - tStart) / (nTime);
//     double time = tStart;
//     HydroForest::RungeKutta3<double> timeIntegrator;
//     for (auto n = 0; n <= nTime; n++) {
//         double dt = timeIntegrator.getMaxTimeStep(2.0, mesh[0].dx());
//         if (time + dt > tEnd) {
//             dt = tEnd - time;
//         }
//         time += dt;
//         app.log("Timestep = %i, Time = %f, dt = %f", n, time, dt);

//         // Construct global solution
//         HydroForest::Vector<double> q_n(nDOFs, 0);
//         for (auto e = 0; e < mesh.size(); e++) {
//             HydroForest::Element1D<double>& element = mesh[e];
//             for (auto i = 0; i < element.size(); i++) {
//                 int I = ID(i,e);
//                 q_n[I] = element.solution()[i];
//             }
//         }

//         // Iterate over elements
//         // HydroForest::Vector<double> R_i;
//         for (auto e = 0; e < mesh.size(); e++) {
//             HydroForest::Element1D<double>& element = mesh[e];

//             // Compute element RHS vector
//             // HydroForest::Matrix<double>& Dtilde_ij = mesh[e].matrices()["weak-derivative-matrix"];
//             HydroForest::Matrix<double>& D_ij = mesh[e].matrices()["strong-derivative-matrix"];
//             HydroForest::Vector<double> f_j = 2.0*mesh[e].solution();
//             element.vectors()["volume-integral"] = -(D_ij * f_j);
//             // R_i = -(D_ij * f_j);

//             // Iterate over nodal points
//             // for (auto i = 0; i < mesh[e].size(); i++) {

//             //     // Update R_I
//             //     int I = ID(i,e);
//             //     R_I[I] += R_i[i];

//             // }

//             // Apply BC


//             // HydroForest::Vector<double> R_j = HydroForest::solve(element.matrices()["mass-matrix"], R_i);

//             // // Update each element solution
//             // // HydroForest::Vector<double> q_next = timeIntegrator.update(time, dt, element.solution(), R_i);
//             // HydroForest::Vector<double> q_next = element.solution() + dt*R_j;
            
//             // element.solution() = q_next;

//         }

//         // Construct global RHS vector
//         HydroForest::Vector<double> R_I = HydroForest::directStiffnessSummationVector(mesh.elements(), ID, "volume-integral");
//         HydroForest::Vector<double> R_J = solve(M_IJ, R_I);
//         // HydroForest::Vector<double> q_update = q_n + dt*R_J;
//         HydroForest::Vector<double> q_update = timeIntegrator.update(time, dt, q_n, R_J);

//         HydroForest::Matrix<double> M_IJ_scaled = (1/(mesh.elements()[0].dx()/6)) * M_IJ;
//         std::cout << "M_IJ = " << M_IJ << std::endl;
//         std::cout << "M_IJ scaled = " << M_IJ_scaled << std::endl;
//         std::cout << "D_IJ = " << mesh[0].matrices()["strong-derivative-matrix"] << std::endl;
//         std::cout << "q_n = " << q_n << std::endl;
//         std::cout << "q_update = " << q_update << std::endl;

//         // Set global solution into mesh
//         for (auto e = 0; e < mesh.size(); e++) {
//             HydroForest::Element1D<double>& element = mesh[e];
//             for (auto i = 0; i < element.size(); i++) {
//                 HydroForest::Vector<int> elementMap = ID.getCol(e);
//                 element.solution()[i] = q_update[elementMap[i]];
//             }
//         }

//         // HydroForest::Vector<double> q_I = HydroForest::directStiffnessSummation(mesh.elements(), ID, )

//         // // Iterate over faces
//         // int nFaces = 2;
//         // for (auto s = 0; s < nFaces; s++) {

//         //     // Apply boundary conditions

//         // }

//         // Construct global R_I

//         // Construct gridpoint solution
//         // HydroForest::Vector<double> R_I_gridpoint = HydroForest::solve(M_IJ, R_I);

//         // 

//         // Update solution via time integration
        
//         // timeIntegrator.update(time, dt, R_I_gridpoint, [&](HydroForest::Vector<double> Rq_n){
//         //     return R_I_gridpoint;
//         // });

//         // Post process and plot
//         mesh.plot("-or");
//         plt::title("Time = " + std::to_string(time));
//         plt::xlim(-1.0, 1.0);
//         plt::ylim(-0.2, 1.2);
//         plt::show();

//         if (time >= tEnd) {
//             break;
//         }

//     }

// }

void runCG(HydroForest::ElementMesh1D<double>& mesh) {

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

    // // Create global grid
    // HydroForest::Vector<double> x(nDOFs);
    // for (auto e = 0; e < nElements; e++) {
    //     HydroForest::Element1D<double>& element = mesh[e];
    //     for (auto i = 0; i < element.size(); i++) {
    //         int I = ID(i,e);
    //         x[I] = element.transformLocal2Global(element.grid()->operator[](i));
    //     }
    // }

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
    Dhat *= 2.0;
    HydroForest::Matrix<double> nDhat = -Dhat;

    std::cout << "--=== MATRICES ===--" << std::endl;
    std::cout << "nElements = " << nElements << "  N = " << N << std::endl;
    std::cout << "M_ij = " << M_ij << std::endl;
    std::cout << "D_ij = " << D_ij << std::endl;
    std::cout << "ID = " << ID << std::endl;
    std::cout << "M_IJ = " << M_IJ << std::endl;
    std::cout << "D_IJ = " << D_IJ << std::endl;
    std::cout << "Dhat = " << Dhat << std::endl;
    std::cout << "-Dhat = " << nDhat << std::endl;

    // Iterate over time
    double tStart = 0;
    double tEnd = 1;
    int nTime = 100000;
    double time = tStart;
    HydroForest::RungeKutta3<double> timeIntegrator;
    std::vector<double> plotTimes = {.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0};
    std::size_t nPlot = 0;
    bool doPlot = false;
    for (auto n = 0; n <= nTime; n++) {
        double dt = 1e-3*timeIntegrator.getMaxTimeStep(2.0, mesh[0].dx());
        double CFL = 2.0*dt/mesh[0].dx();
        // double dt = 0.01;
        if (time + dt > tEnd) {
            dt = tEnd - time;
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
        // std::cout << "q_n = " << q_n << std::endl;

        // HydroForest::Vector<double> f_n = 2.0*q_n;
        HydroForest::Vector<double> q_update = 1.0*timeIntegrator.update(time, dt, q_n, -Dhat);
        // for (auto i = 0; i < q_update.size(); i++) {
        //     q_update[i] = -1.0*q_update[i];
        // }
        // HydroForest::Vector<double> r_n = Dhat*q_n;
        // HydroForest::Vector<double> q_update = q_n - 2.0*dt*r_n;
        // std::cout << "q_n = " << q_n << std::endl;
        // std::cout << "q_update = " << q_update << std::endl; 

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
            // plt::plot(x.data(), q_update.data(), "-or");
            plt::title("Time = " + std::to_string(time));
            plt::xlim(-1.0, 1.0);
            plt::ylim(-0.2, 1.2);
            plt::show();

            nPlot++;
            doPlot = false;
        }

        if (time >= tEnd) {
            break;
        }
    }

}

void runDG(HydroForest::ElementMesh1D<double>& mesh) {

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

    // Iterate over time
    double tStart = 0;
    double tEnd = 1;
    int nTime = 100000;
    double time = tStart;
    HydroForest::RungeKutta3<double> timeIntegrator;
    std::vector<double> plotTimes = {.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0};
    std::size_t nPlot = 0;
    bool doPlot = false;
    for (auto n = 0; n <= nTime; n++) {
        double dt = 1e-1*timeIntegrator.getMaxTimeStep(2.0, mesh[0].dx());
        double CFL = 2.0*dt/mesh[0].dx();
        if (time + dt > tEnd) {
            dt = tEnd - time;
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

        if (time >= tEnd) {
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
    std::cout << mesh << std::endl;

    // Plot initial condition
    mesh.plot("-or");
    plt::title("Initial Condition");
    plt::xlim(-1.0, 1.0);
    plt::ylim(-0.2, 1.2);
    plt::show();

    if (std::get<std::string>(options["scheme"]) == "CG") {
        runCG(mesh);
    }
    else if (std::get<std::string>(options["scheme"]) == "DG") {
        runDG(mesh);
    }
    else {
        throw std::invalid_argument("Invalid Scheme; Options are `CG` or `DG`");
    }

    return EXIT_SUCCESS;
}