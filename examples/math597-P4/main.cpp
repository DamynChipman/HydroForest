#include <cmath>
#include <iostream>

#include <p4est.h>
#include <p4est_connectivity.h>
#include <p4est_vtk.h>

#include <P4est.hpp>
#include <Grid2D.hpp>
#include <DGFiniteElements.hpp>
#include <VTKBuilder.hpp>
#include <Element2D.hpp>
#include <Mesh2D.hpp>

// static auto f = [](double x, double y) { return sin(M_PI*x) * cos(M_PI*y); };
static auto f = [](double x, double y) { return x*x + y*x + y*y; };

void plotBasisFunctions(int xOrder, int yOrder, int xPlot, int yPlot) {
    HydroForest::LobattoTensorProductGrid2D<double> grid(xOrder, yOrder);
    HydroForest::LobattoTensorProductGrid2D<double> quadratureGrid(xOrder+1, yOrder+1);
    HydroForest::LobattoTensorProductGrid2D<double> plotGrid(xPlot, yPlot);
    HydroForest::Matrix<double> L_I = grid.basisMatrix(plotGrid);

    std::vector<HydroForest::Vector<double>> psi_I(grid.size());
    for (auto i = 0; i < grid.xSize(); i++) {
        for (auto j = 0; j < grid.ySize(); j++) {
            int I = grid.ID(i,j);
            psi_I[I] = L_I.getRow(I);
            psi_I[I].setName(std::to_string(I));

        }
    }

    HydroForest::RectilinearGridVTK rectilinearGridVTK;
    rectilinearGridVTK.buildMesh(plotGrid, &plotGrid.xGrid(), &plotGrid.yGrid(), nullptr);
    for (auto& psi : psi_I) rectilinearGridVTK.addPointData(psi);
    rectilinearGridVTK.toVTK("basis_functions.vtr");
}

int main(int argc, char** argv) {

    HydroForest::HydroForestApp app(&argc, &argv);

    // Create p4est mesh
    int nElementsPerSide = 2;
    int order = 1;
    double xLower = -2;
    double xUpper = 2;
    double yLower = -2;
    double yUpper = 2;
    HydroForest::ElementMesh2D<double> mesh(xLower, xUpper, yLower, yUpper, nElementsPerSide, order);

    // mesh.iterate([&](HydroForest::QuadElement2D<double>& element){
    //     std::cout << "element = " << element.ID() << std::endl;
    //     for (auto l = 0; l < 4; l++) {
    //         std::cout << "  l = " << l << " : neighbor  = " << element.neighbors()[l]->ID() << std::endl;
    //     }
    // });

    // Physical quantities
    double x0 = -0.5;
    double y0 = 0.0;
    double sigma = 1.0/8.0;
    double speed = 1.0;

    // Set initial condition on mesh
    mesh.setVariables({"q", "ux", "uy", "fx", "fy"});
    mesh.setSolution("q", [&](double x, double y){
        return exp(-(pow(x - x0, 2) + pow(y - y0, 2)) / (2.0*pow(sigma, 2)));
    });
    mesh.setSolution("ux", [&](double x, double y){
        return y;
    });
    mesh.setSolution("uy", [&](double x, double y){
        return -x;
    });
    mesh.iterate([&](HydroForest::QuadElement2D<double>& element){
        HydroForest::Vector<double>& q = element.vector("q");
        HydroForest::Vector<double>& ux = element.vector("ux");
        HydroForest::Vector<double>& uy = element.vector("uy");
        HydroForest::Vector<double>& fx = element.vector("fx");
        HydroForest::Vector<double>& fy = element.vector("fy");
        
        for (auto i = 0; i < fx.size(); i++) fx[i] = q[i] * ux[i];
        for (auto i = 0; i < fy.size(); i++) fy[i] = q[i] * uy[i];
    });
    mesh.toVTK("initial_condition", {"q", "ux", "uy", "fx", "fy"});

    HydroForest::QuadElement2D<double>& element0 = mesh.element(3);
    HydroForest::Vector<double> R_volume = element0.computeVolumeIntegralInexact();
    HydroForest::Vector<double> R_flux = element0.computeFluxIntegralInexact(0);
    std::cout << "R_volume = " << R_volume << std::endl;
    std::cout << "R_flux = " << R_flux << std::endl;
    // HydroForest::Vector<double> R_flux = mesh.element(0).computeFluxIntegralInexact(0);

    // Set matrices on mesh
    // mesh.setMassMatrices();
    // mesh.setDerivativeMatrices();

    // Test element
    // std::vector<HydroForest::SpaceVector2D<double>> points = {
    //     HydroForest::SpaceVector2D<double>(-1, -1),
    //     HydroForest::SpaceVector2D<double>(1, -1),
    //     HydroForest::SpaceVector2D<double>(1, 1),
    //     HydroForest::SpaceVector2D<double>(-1, 1)
    // };
    // HydroForest::QuadElement2D<double> quad(points, order, order, -1);

    // quad.matrixMap()["mass"] = HydroForest::DGMassMatrix2D<double>(quad);
    // // quad.matrixMap()["derivative-x"] = HydroForest::DGDerivativeMatrix2D<double>(quad, 0);
    // // quad.matrixMap()["derivative-y"] = HydroForest::DGDerivativeMatrix2D<double>(quad, 1);
    // // quad.matrixMap()["flux-x-0"] = HydroForest::DGFluxMatrix2D<double>(quad, 0, 0);
    // // quad.matrixMap()["flux-y-0"] = HydroForest::DGFluxMatrix2D<double>(quad, 1, 0);
    // // quad.matrixMap()["flux-x-1"] = HydroForest::DGFluxMatrix2D<double>(quad, 0, 1);
    // // quad.matrixMap()["flux-y-1"] = HydroForest::DGFluxMatrix2D<double>(quad, 1, 1);
    // // quad.matrixMap()["flux-x-2"] = HydroForest::DGFluxMatrix2D<double>(quad, 0, 2);
    // // quad.matrixMap()["flux-y-2"] = HydroForest::DGFluxMatrix2D<double>(quad, 1, 2);
    // // quad.matrixMap()["flux-x-3"] = HydroForest::DGFluxMatrix2D<double>(quad, 0, 3);
    // // quad.matrixMap()["flux-y-3"] = HydroForest::DGFluxMatrix2D<double>(quad, 1, 3);

    // if (order == 1 || order == 2) {

    //     std::cout << "jacobian = " << quad.vector("jacobian") << std::endl;

    //     for (auto l = 0; l < 4; l++) {
    //         std::cout << "l = " << l << std::endl;
    //         std::cout << "normals = " << quad.normals()[l] << std::endl;
    //     }

    //     double dx = 2.0;
    //     double dy = 2.0;
    //     quad.matrix("mass") *= (36.0 / (dx*dy));
    //     // quad.matrix("derivative-x") *= (12.0 / dy);
    //     // quad.matrix("derivative-y") *= (12.0 / dx);
    //     // quad.matrix("flux-x-0") *= (6.0 / dx);
    //     // quad.matrix("flux-y-0") *= (6.0 / dx);
    //     // quad.matrix("flux-x-1") *= (6.0 / dy);
    //     // quad.matrix("flux-y-1") *= (6.0 / dy);
    //     // quad.matrix("flux-x-2") *= (6.0 / dx);
    //     // quad.matrix("flux-y-2") *= (6.0 / dx);
    //     // quad.matrix("flux-x-3") *= (6.0 / dy);
    //     // quad.matrix("flux-y-3") *= (6.0 / dy);

    //     std::cout << "M = " << quad.matrix("mass") << std::endl;
    //     // std::cout << "Dx = " << quad.matrix("derivative-x") << std::endl;
    //     // std::cout << "Dy = " << quad.matrix("derivative-y") << std::endl;
    //     // std::cout << "flux-x-0 = " << quad.matrix("flux-x-0") << std::endl;
    //     // std::cout << "flux-y-0 = " << quad.matrix("flux-y-0") << std::endl;
    //     // std::cout << "flux-x-1 = " << quad.matrix("flux-x-1") << std::endl;
    //     // std::cout << "flux-y-1 = " << quad.matrix("flux-y-1") << std::endl;
    //     // std::cout << "flux-x-2 = " << quad.matrix("flux-x-2") << std::endl;
    //     // std::cout << "flux-y-2 = " << quad.matrix("flux-y-2") << std::endl;
    //     // std::cout << "flux-x-3 = " << quad.matrix("flux-x-3") << std::endl;
    //     // std::cout << "flux-y-3 = " << quad.matrix("flux-y-3") << std::endl;
    // }
    
    // M *= (9.0);
    // std::cout << "M = " << M << std::endl;
    // std::cout << "D = " << D << std::endl;

    return 0;
}