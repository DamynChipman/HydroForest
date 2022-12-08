#ifndef MESH_2D_HPP_
#define MESH_2D_HPP_

#include "P4est.hpp"
#include "Element2D.hpp"
#include "VTKBuilder.hpp"

namespace HydroForest {

template<typename FloatingDataType>
class ElementMesh2D {

protected:

    FloatingDataType xLower_;
    FloatingDataType xUpper_;
    FloatingDataType yLower_;
    FloatingDataType yUpper_;
    int nElements_;
    int order_;
    p4est_t* p4est_;
    std::vector<QuadElement2D<FloatingDataType>> elements_;

public:

    ElementMesh2D(FloatingDataType xLower, FloatingDataType xUpper, FloatingDataType yLower, FloatingDataType yUpper, int nElementsPerSide, int order) :
        xLower_(xLower), xUpper_(xUpper), yLower_(yLower), yUpper_(yUpper), nElements_(nElementsPerSide*nElementsPerSide), order_(order) {

        // Create p4est mesh
        p4est_connectivity_t* initialConnectivity = p4est::p4est_connectivity_new_rectangular_domain(xLower, xUpper, yLower, yUpper);
        p4est_connectivity_t* refinedConnectivity = p4est_connectivity_refine(initialConnectivity, nElementsPerSide);
        p4est_ = p4est_new(MPI_COMM_WORLD, refinedConnectivity, sizeof(ElementMesh2D), p4est_element_init_, this);

        // Write initial mesh to file
        p4est_vtk_write_file(p4est_, NULL, "initial_mesh");

    }

    FloatingDataType xLower() { return xLower_; }
    FloatingDataType xUpper() { return xUpper_; }
    FloatingDataType yLower() { return yLower_; }
    FloatingDataType yUpper() { return yUpper_; }
    int size() { return nElements_; }
    int order() { return order_; }
    std::vector<QuadElement2D<FloatingDataType>>& elements() { return elements_; }

    void setVariables(std::vector<std::string> variableNames) {
        for (auto e = 0; e < nElements_; e++) {
            QuadElement2D<FloatingDataType>& element = elements_[e];
            for (auto s = 0; s < variableNames.size(); s++) {
                element.vectorMap()[variableNames[s]] = Vector<FloatingDataType>(element.referenceGrid().size(), 0);
            }
        }
    }

    void setSolution(std::string variableName, std::function<FloatingDataType(FloatingDataType, FloatingDataType)> function) {

        for (auto e = 0; e < nElements_; e++) {
            QuadElement2D<FloatingDataType>& element = elements_[e];
            LobattoTensorProductGrid2D<FloatingDataType>& grid = element.referenceGrid();
            for (auto i = 0; i < grid.size(); i++) {
                SpaceVector2D<FloatingDataType> xi = grid(i);
                SpaceVector2D<FloatingDataType> x = element.mapReference2Physical(xi);
                element.vectorMap()[variableName][i] = function(x.x(), x.y());
            }
        }

    }

    void toVTK(std::string filename, std::vector<std::string> variableNames) {

        for (auto e = 0; e < nElements_; e++) {
            QuadElement2D<FloatingDataType>& element = elements_[e];
            LobattoTensorProductGrid2D<FloatingDataType>& grid = element.physicalGrid();

            RectilinearGridVTK gridVTK;
            gridVTK.buildMesh(grid, &grid.xGrid(), &grid.yGrid(), nullptr);
            for (auto s = 0; s < variableNames.size(); s++) {
                element.vector(variableNames[s]).setName(variableNames[s]);
                gridVTK.addPointData(element.vector(variableNames[s]));
            }
            gridVTK.toVTK(filename + "_" + std::to_string(e) + ".vtr");

        }

    }

private:

    static void p4est_element_init_(p4est_t* p4est, p4est_topidx_t which_tree, p4est_quadrant_t* quadrant) {
        
        // Get mesh and elements
        ElementMesh2D<FloatingDataType>& mesh = *((ElementMesh2D<FloatingDataType>*) p4est->user_pointer);
        std::vector<QuadElement2D<FloatingDataType>>& elements = mesh.elements();

        // Get quadrant points
        p4est_qcoord_t quad_length = P4EST_QUADRANT_LEN(quadrant->level);
        double xyz0[3];
        double xyz1[3];
        double xyz2[3];
        double xyz3[3];
        p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x, quadrant->y, xyz0);
        p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x + quad_length, quadrant->y, xyz1);
        p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x + quad_length, quadrant->y + quad_length, xyz2);
        p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x, quadrant->y + quad_length, xyz3);
        
        // Create vector of points
        std::vector<SpaceVector2D<FloatingDataType>> points = {
            SpaceVector2D<FloatingDataType>(xyz0[0], xyz0[1]),
            SpaceVector2D<FloatingDataType>(xyz1[0], xyz1[1]),
            SpaceVector2D<FloatingDataType>(xyz2[0], xyz2[1]),
            SpaceVector2D<FloatingDataType>(xyz3[0], xyz3[1])
        };

        // for (auto i = 0; i < 4; i++) {
        //     std::cout << "point " << i << ": " << "[" << points[i].x() << ", " << points[i].y() << "]" << std::endl;
        // }

        // Create element and add to elements
        elements.push_back(QuadElement2D<FloatingDataType>(points, mesh.order(), mesh.order()));

    }

};

} // NAMESPACE : HydroForest

#endif // MESH_2D_HPP_