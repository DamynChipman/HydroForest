#include <cmath>
#include <iostream>

#include <HydroForestApp.hpp>
#include <Mesh1D.hpp>
#include <CGFiniteElements.hpp>
#include <DGFiniteElements.hpp>

int PLOT_ID = 2;

double f_RHS(double x) {
    return (1.0 - pow(M_PI,2))*sin(M_PI*x);
}

double q_exact(double x) {
    return sin(M_PI*x);
}

double computeError(HydroForest::Vector<double>& qNumerical, HydroForest::Vector<double>& qExact) {

    if (qNumerical.size() != qExact.size()) {
        throw std::invalid_argument("Sizes of `qNumerical` and `qExact` are not the same");
    }

    double numerator = 0.0;
    double denominator = 0.0;
    for (auto k = 0; k < qNumerical.size(); k++) {
        numerator += pow(qNumerical[k] - qExact[k], 2.0);
        denominator += pow(qExact[k], 2.0);
    }

    return pow(numerator / denominator, 0.5);

}

double solveCG(HydroForest::ElementMesh1D<double>& mesh, bool plotFlag) {

    // App and options
    HydroForest::HydroForestApp& app = HydroForest::HydroForestApp::getInstance();
    HydroForest::Options& options = app.getOptions();

    // Get options
    std::string scheme = std::get<std::string>(options["scheme"]);
    std::string integration = std::get<std::string>(options["integration"]);

    // Get values from the mesh
    int nOrder = mesh.order();
    int nElements = mesh.size();
    int nDOFs = nOrder*nElements + 1;
    std::vector<HydroForest::Element1D<double>>& elements = mesh.elements();
    HydroForest::LagrangePolynomial& poly = elements[0].polynomial();
    HydroForest::Grid1DBase<double>& grid = *elements[0].grid();
    HydroForest::Grid1DBase<double>& qGrid = *elements[0].quadratureGrid();

    // Create element matrices
    HydroForest::CGMassMatrix<double> M_ij(poly, grid, qGrid);
    HydroForest::CGLaplacianMatrix<double> L_ij(poly, grid, qGrid);

    // if (nOrder == 1) {
    //     std::cout << "M_ij = " << M_ij << std::endl;
    //     std::cout << "L_ij = " << L_ij << std::endl;
    // }

    // Create ID matrix
    HydroForest::CGIDMatrix ID(nOrder, nElements);

    // Create global matrices
    HydroForest::CGDirectStiffnessSummationOperator<double> DSS(elements, ID);
    HydroForest::Matrix<double> M_IJ = DSS.operateWithMetricTerm(M_ij);
    HydroForest::Matrix<double> L_IJ = DSS.operateWithInverseMetricTerm(L_ij);
    HydroForest::Matrix<double> ML_IJ = M_IJ - L_IJ;

    // if (nOrder == 1) {
    //     std::cout << "M_IJ = " << M_IJ << std::endl;
    //     std::cout << "L_IJ = " << L_IJ << std::endl;
    //     std::cout << "M_IJ - L_IJ = " << ML_IJ << std::endl;
    // }

    // Create global RHS function
    HydroForest::Vector<double> f_J(nDOFs, 0);
    for (auto e = 0; e < mesh.size(); e++) {
        HydroForest::Element1D<double>& element = mesh[e];
        for (auto i = 0; i < element.size(); i++) {
            double xi = element.grid()->operator[](i);
            double x = element.transformLocal2Global(xi);
            int I = ID(i,e);
            f_J[I] = f_RHS(x);
        }
    }

    // Apply BC
    // Modify LHS matrix for left side Dirichlet BC
    HydroForest::Vector<double> zeroRow(ML_IJ.nCols(), 0);
    ML_IJ.setRow(0, zeroRow);
    ML_IJ(0,0) = 1.0;

    // if (nOrder == 1) {
    //     std::cout << "M_IJ - L_IJ = " << ML_IJ << std::endl;
    // }

    // Neumann boundary vector
    HydroForest::Vector<double> B_I(nDOFs, 0);
    B_I[nDOFs-1] = -M_PI;

    // Compute intermediate vectors
    HydroForest::Vector<double> R_I = M_IJ*f_J;
    HydroForest::Vector<double> RHS_I = R_I - B_I;

    // Solve system
    HydroForest::Vector<double> q_J = HydroForest::solve(ML_IJ, RHS_I);

    // if (nOrder == 1) {
    //     std::cout << "q_J = " << q_J << std::endl;
    // }

    // Set solution into mesh for plotting
    for (auto e = 0; e < mesh.size(); e++) {
        HydroForest::Element1D<double>& element = mesh[e];
        for (auto i = 0; i < element.size(); i++) {
            HydroForest::Vector<int> elementMap = ID.getCol(e);
            element.solution()[i] = q_J[elementMap[i]];
        }
    }

    // Create global exact solution
    HydroForest::Vector<double> x_J(nDOFs, 0);
    HydroForest::Vector<double> q_exact_J(nDOFs, 0);
    for (auto e = 0; e < mesh.size(); e++) {
        HydroForest::Element1D<double>& element = mesh[e];
        for (auto i = 0; i < element.size(); i++) {
            double xi = element.grid()->operator[](i);
            double x = element.transformLocal2Global(xi);
            int I = ID(i,e);
            q_exact_J[I] = q_exact(x);
            x_J[I] = x;
        }
    }

    // Compute error
    double error = computeError(q_J, q_exact_J);
    app.log("nElements = %i,  nOrder = %i,  nDOFs = %i,  L2-error = %e", nElements, nOrder, nDOFs, error);

    // Plot
    if (plotFlag) {
        plt::figure(PLOT_ID++);
        mesh.plot("-or");
        plt::named_plot("Exact Solution", x_J.data(), q_exact_J.data(), "-b");
        plt::title("Helmholtz Equation - " + scheme + "-" + integration + ":\nN = " + std::to_string(nOrder) + ",  nElements = " + std::to_string(nElements) + ",  nDOFs = " + std::to_string(nDOFs));
        plt::legend();
        plt::save("plot_solution_" + scheme + "_" + integration + "_p" + std::to_string((int) nOrder) + "_n" + std::to_string((int) nElements) + ".png");
        // plt::show();
    }

    return error;

}

double solveDG(HydroForest::ElementMesh1D<double>& mesh, bool plotFlag) {

    // App and options
    HydroForest::HydroForestApp& app = HydroForest::HydroForestApp::getInstance();
    HydroForest::Options& options = app.getOptions();

    // Get options
    std::string scheme = std::get<std::string>(options["scheme"]);
    std::string integration = std::get<std::string>(options["integration"]);

    // Get values from the mesh
    int nOrder = mesh.order();
    int nElements = mesh.size();
    int nDOFs = nElements*(nOrder + 1);
    std::vector<HydroForest::Element1D<double>>& elements = mesh.elements();
    HydroForest::LagrangePolynomial& poly = elements[0].polynomial();
    HydroForest::Grid1DBase<double>& grid = *elements[0].grid();
    HydroForest::Grid1DBase<double>& qGrid = *elements[0].quadratureGrid();

    // Create element matrices
    HydroForest::DGMassMatrix<double> M_ij(poly, grid, qGrid);
    HydroForest::DGDerivativeMatrix<double> Dtilde_ij(poly, grid, qGrid);
    HydroForest::DGFluxMatrix<double> F_ij(nOrder);

    // if (nOrder == 1) {
    //     std::cout << "M_ij = " << M_ij << std::endl;
    //     std::cout << "Dtilde_ij = " << Dtilde_ij << std::endl;
    // }

    // Create ID matrix
    HydroForest::DGIDMatrix ID(nOrder, nElements);

    // Create global matrices
    HydroForest::DGDirectStiffnessSummationOperator<double> DSS(elements, ID);
    HydroForest::Matrix<double> M_IJ = DSS.operateWithMetricTerm(M_ij);
    HydroForest::Matrix<double> Dtilde_IJ = DSS.operate(Dtilde_ij);

    // Create global flux matrices
    HydroForest::DGGlobalCenteredFluxMatrix<double> F_IJ_q(elements, ID, HydroForest::BoundaryConditionType::Dirichlet, HydroForest::BoundaryConditionType::Neumann);
    HydroForest::DGGlobalCenteredFluxMatrix<double> F_IJ_Q(elements, ID, HydroForest::BoundaryConditionType::Dirichlet, HydroForest::BoundaryConditionType::Neumann);

    // Include the "action" of B_Q and B_q to F
    F_IJ_q(0,0) = -0.5;
    F_IJ_q(nDOFs-1,nDOFs-1) = 1;
    F_IJ_Q(0,0) = -1;
    F_IJ_Q(nDOFs-1,nDOFs-1) = 0.5;

    // Create intermediate matrices
    HydroForest::Matrix<double> Dhat_IJ_q = F_IJ_q - Dtilde_IJ;
    HydroForest::Matrix<double> Dhat_IJ_Q = F_IJ_Q - Dtilde_IJ;

    // Create Helmholtz matrix
    HydroForest::Matrix<double> H_IJ = HydroForest::solve(M_IJ, Dhat_IJ_q);
    H_IJ = Dhat_IJ_Q * H_IJ;
    H_IJ = H_IJ + M_IJ;

    // Create boundary condition vectors
    HydroForest::Vector<double> B_I_q(nDOFs, 0);
    B_I_q[0] = 0.0; // -0.5*g(xL)
    HydroForest::Vector<double> B_I_Q(nDOFs, 0);
    B_I_Q[nDOFs-1] = -0.5*M_PI; // 0.5*h(xR)
    HydroForest::Vector<double> Bprime_I = HydroForest::solve(M_IJ, B_I_q);
    Bprime_I = Dhat_IJ_Q * Bprime_I;
    Bprime_I = B_I_Q - Bprime_I;

    // Create global RHS function
    HydroForest::Vector<double> f_J(nDOFs, 0);
    for (auto e = 0; e < mesh.size(); e++) {
        HydroForest::Element1D<double>& element = mesh[e];
        for (auto i = 0; i < element.size(); i++) {
            double xi = element.grid()->operator[](i);
            double x = element.transformLocal2Global(xi);
            int I = ID(i,e);
            f_J[I] = f_RHS(x);
        }
    }

    // Create RHS vector
    HydroForest::Vector<double> RHS_I = M_IJ*f_J;
    RHS_I = RHS_I - Bprime_I;

    // Solve system
    HydroForest::Vector<double> q_J = HydroForest::solve(H_IJ, RHS_I);

    // if (nOrder <= 2) {
    //     std::cout << "M_IJ = " << M_IJ << std::endl;
    //     std::cout << "Dtilde_IJ = " << Dtilde_IJ << std::endl;
    //     std::cout << "F_IJ_q = " << F_IJ_q << std::endl;
    //     std::cout << "F_IJ_Q = " << F_IJ_Q << std::endl;
    //     std::cout << "Dhat_IJ_q = " << Dhat_IJ_q << std::endl;
    //     std::cout << "Dhat_IJ_Q = " << Dhat_IJ_Q << std::endl;
    //     std::cout << "H_IJ = " << H_IJ << std::endl;
    //     std::cout << "Bprime_I = " << Bprime_I << std::endl;
    //     std::cout << "RHS_I = " << RHS_I << std::endl;
    //     std::cout << "q_J = " << q_J << std::endl;
    // }


    // Set solution into mesh for plotting
    for (auto e = 0; e < mesh.size(); e++) {
        HydroForest::Element1D<double>& element = mesh[e];
        for (auto i = 0; i < element.size(); i++) {
            HydroForest::Vector<int> elementMap = ID.getCol(e);
            element.solution()[i] = q_J[elementMap[i]];
        }
    }

    // Create global exact solution
    HydroForest::Vector<double> x_J(nDOFs, 0);
    HydroForest::Vector<double> q_exact_J(nDOFs, 0);
    for (auto e = 0; e < mesh.size(); e++) {
        HydroForest::Element1D<double>& element = mesh[e];
        for (auto i = 0; i < element.size(); i++) {
            double xi = element.grid()->operator[](i);
            double x = element.transformLocal2Global(xi);
            int I = ID(i,e);
            q_exact_J[I] = q_exact(x);
            x_J[I] = x;
        }
    }

    // Compute error
    double error = computeError(q_J, q_exact_J);
    app.log("nElements = %i,  nOrder = %i,  nDOFs = %i,  L2-error = %e", nElements, nOrder, nDOFs, error);

    // Plot
    if (plotFlag) {
        plt::figure(PLOT_ID++);
        mesh.plot("-or");
        plt::named_plot("Exact Solution", x_J.data(), q_exact_J.data(), "-b");
        plt::title("Helmholtz Equation - " + scheme + "-" + integration + ":\nN = " + std::to_string(nOrder) + ",  nElements = " + std::to_string(nElements) + ",  nDOFs = " + std::to_string(nDOFs));
        plt::legend();
        plt::save("plot_solution_" + scheme + "_" + integration + "_p" + std::to_string((int) nOrder) + "_n" + std::to_string((int) nElements) + ".png");
        // plt::show();
    }

    return error;

}

int main(int argc, char** argv) {

    HydroForest::HydroForestApp app(&argc, &argv);
    HydroForest::Options& options = app.getOptions();

    app.log("Hello from P3!");
    app.log("Solving Helmholtz Equation:");
    app.log("  PDE:  d2q(x)/dx2 + q(x) = (1 - M_PI^2) sin(M_PI*x)");
    app.log("  BC:   q(-1) = 0");
    app.log("  BC:   dq/dx|x=1 = -M_PI");

    // Setup options to run through
    std::vector<std::string> schemeVector = {"CG", "DG"};
    std::vector<std::string> integrationVector = {"exact", "inexact"};
    std::vector<double> nOrderVector = {1, 2, 4, 8, 16, 32};
    std::vector<double> nElementVector = {4, 8, 16, 32, 64, 128};

    int mpiSize;
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    if (mpiSize != 4) {
        std::cerr << "Must be run with 4 ranks." << std::endl;
        return EXIT_FAILURE;
    }

    if (myRank == 0) {
        options["scheme"] = "CG";
        options["integration"] = "exact";
    }
    else if (myRank == 1) {
        options["scheme"] = "CG";
        options["integration"] = "inexact";
    }
    else if (myRank == 2) {
        options["scheme"] = "DG";
        options["integration"] = "exact";
    }
    else if (myRank == 3) {
        options["scheme"] = "DG";
        options["integration"] = "inexact";
    }
    std::string scheme = std::get<std::string>(options["scheme"]);
    std::string integration = std::get<std::string>(options["integration"]);
    int plotFinalID = 2;

    for (auto& nOrder : nOrderVector) {
        std::vector<double> errorVector;
        std::vector<double> nDOFsVector;
        for (auto& nElement : nElementVector) {

            // Create the mesh
            double xLower = -1.0;
            double xUpper = 1.0;
            HydroForest::ElementMesh1D<double> mesh(xLower, xUpper, nElement, nOrder);

            // Solve with CG or DG methods and get error
            double error;
            int nDOFs;
            bool plotFlag = true;
            if (std::get<std::string>(options["scheme"]) == "CG") {
                error = solveCG(mesh, plotFlag);
                nDOFs = nOrder*nElement + 1;
            }
            else if (std::get<std::string>(options["scheme"]) == "DG") {
                error = solveDG(mesh, plotFlag);
                nDOFs = nElement*(nOrder + 1);
            }
            else {
                throw std::invalid_argument("Invalid scheme; Options are `CG` or `DG`");
            }

            // Store error and nDOFs
            errorVector.push_back(error);
            nDOFsVector.push_back(nDOFs);

        }

        plt::figure(1);
        plt::named_loglog("p = " + std::to_string((int) nOrder), nDOFsVector, errorVector, "-*");
    }

    // Make the plot pretty
    std::vector<int> xTicks = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
    std::vector<std::string> xTickLabels;
    for (auto& t : xTicks) xTickLabels.push_back(std::to_string(t));
    plt::title(scheme + "-" + integration);
    plt::xlabel("Total Degrees of Freedom");
    plt::ylabel("${L_2}$ Error");
    plt::xticks(xTicks, xTickLabels);
    plt::legend();
    plt::grid(true);
    plt::save("plot_convergance_" + scheme + "_" + integration + ".png");

    return EXIT_SUCCESS;
}