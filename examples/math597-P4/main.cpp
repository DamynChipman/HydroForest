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

static auto f = [](double x, double y) { return sin(M_PI*x) * cos(M_PI*y); };

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

    // Set initial condition on mesh
    mesh.setVariables({"q"});
    mesh.setSolution("q", [](double x, double y){
        return x + y;
    });
    mesh.toVTK("initial_condition", {"q"});
    
    int xOrder = 1;
    int yOrder = 1;
    int xPlot = 50;
    int yPlot = 50;
    HydroForest::LobattoTensorProductGrid2D<double> grid(xOrder, yOrder);
    HydroForest::LobattoTensorProductGrid2D<double> quadratureGrid(xOrder+1, yOrder+1);
    HydroForest::LobattoTensorProductGrid2D<double> plotGrid(xPlot, yPlot);
    HydroForest::Matrix<double> L_I = grid.basisMatrix(plotGrid);
    // std::cout << "L_I = " << L_I << std::endl;
    // HydroForest::LagrangePolynomial psi_x(grid.xPoints());
    // HydroForest::LagrangePolynomial psi_y(grid.yPoints());
    // HydroForest::Matrix<double> L_x = psi_x(plotGrid.xPoints());
    // HydroForest::Matrix<double> L_y = psi_y(plotGrid.yPoints());

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
    rectilinearGridVTK.toVTK("test.vtr");

    HydroForest::DGMassMatrix2D M(grid, quadratureGrid);
    HydroForest::DGDerivativeMatrix2D D(grid, quadratureGrid);
    // M *= (9.0);
    std::cout << "M = " << M << std::endl;
    std::cout << "D = " << D << std::endl;
    
    std::vector<HydroForest::SpaceVector2D<double>> points = {
        HydroForest::SpaceVector2D<double>(-1, -1),
        HydroForest::SpaceVector2D<double>(1, -1),
        HydroForest::SpaceVector2D<double>(1, 1),
        HydroForest::SpaceVector2D<double>(-1, 1)
    };
    HydroForest::QuadElement2D<double> quad(points, 1, 1);
    std::cout << "jacobian = " << quad.vector("jacobian") << std::endl;

    return 0;
}