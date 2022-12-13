#ifndef MESH_2D_HPP_
#define MESH_2D_HPP_

#include "P4est.hpp"
#include "Element2D.hpp"
#include "VTKBuilder.hpp"
#include "XMLTree.hpp"

namespace HydroForest {

static int ELEMENT_COUNTER = 0;

template<typename FloatingDataType>
using ElementCallbackFunctionType = std::function<void(QuadElement2D<FloatingDataType>&)>;

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
    // std::vector<QuadElement2D<FloatingDataType>> elements_;
    // std::vector<QuadElement2D<FloatingDataType>*> elementsPtrs_;

public:

    ElementMesh2D(FloatingDataType xLower, FloatingDataType xUpper, FloatingDataType yLower, FloatingDataType yUpper, int minLevel, int order) :
        xLower_(xLower), xUpper_(xUpper), yLower_(yLower), yUpper_(yUpper), order_(order) {

        // Create p4est mesh
        // p4est_connectivity_t* initialConnectivity = p4est_connectivity_new_periodic();
        p4est_connectivity_t* initialConnectivity = p4est::p4est_connectivity_new_rectangular_domain(xLower, xUpper, yLower, yUpper);
        // p4est_connectivity_destroy(initialConnectivity);
        // p4est_ = p4est_new(
        //     MPI_COMM_WORLD,
        //     refinedConnectivity,
        //     // initialConnectivity,
        //     sizeof(ElementMesh2D),
        //     [](p4est_t* p4est, p4est_topidx_t which_tree, p4est_quadrant_t* quadrant) {
        
        //         // Get mesh and elements
        //         ElementMesh2D<FloatingDataType>& mesh = *((ElementMesh2D<FloatingDataType>*) p4est->user_pointer);
        //         std::vector<QuadElement2D<FloatingDataType>>& elements = mesh.elements();

        //         // Get quadrant points
        //         p4est_qcoord_t quad_length = P4EST_QUADRANT_LEN(quadrant->level);
        //         double xyz0[3];
        //         double xyz1[3];
        //         double xyz2[3];
        //         double xyz3[3];
        //         p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x, quadrant->y, xyz0);
        //         p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x + quad_length, quadrant->y, xyz1);
        //         p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x + quad_length, quadrant->y + quad_length, xyz2);
        //         p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x, quadrant->y + quad_length, xyz3);
                
        //         // Create vector of points
        //         std::vector<SpaceVector2D<FloatingDataType>> points = {
        //             SpaceVector2D<FloatingDataType>(xyz0[0], xyz0[1]),
        //             SpaceVector2D<FloatingDataType>(xyz1[0], xyz1[1]),
        //             SpaceVector2D<FloatingDataType>(xyz2[0], xyz2[1]),
        //             SpaceVector2D<FloatingDataType>(xyz3[0], xyz3[1])
        //         };

        //         // Create element and add to elements
        //         // QuadElement2D<FloatingDataType>* newElement = new QuadElement2D<FloatingDataType>(points, mesh.order(), mesh.order(), ELEMENT_COUNTER);
        //         elements.emplace_back(points, mesh.order(), mesh.order(), ELEMENT_COUNTER);

        //         // Set user data in p4est quadrant
        //         // quadrant->p.user_data = (void*) newElement;
        //         // quadrant->p.user_data = (void*) &(elements[ELEMENT_COUNTER]);
        //         quadrant->p.user_int = ELEMENT_COUNTER;
        //         // std::cout << "quadrant->p.user_data = " << quadrant->p.user_data << std::endl;

        //         // Update counter
        //         ELEMENT_COUNTER++;

        //         return;

        //     },
        //     this
        // );
        p4est_ = p4est_new_ext(
            MPI_COMM_WORLD,
            // refinedConnectivity,
            initialConnectivity,
            0,
            minLevel,
            1,
            sizeof(ElementMesh2D),
            [](p4est_t* p4est, p4est_topidx_t which_tree, p4est_quadrant_t* quadrant) {
        
                // Get mesh and elements
                ElementMesh2D<FloatingDataType>& mesh = *((ElementMesh2D<FloatingDataType>*) p4est->user_pointer);
                // std::vector<QuadElement2D<FloatingDataType>>& elements = mesh.elements();

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

                // Create element and add to elements
                QuadElement2D<FloatingDataType>* newElement = new QuadElement2D<FloatingDataType>(points, mesh.order(), mesh.order(), ELEMENT_COUNTER);
                // elements.push_back(*newElement);
                // elements.emplace_back(points, mesh.order(), mesh.order(), ELEMENT_COUNTER);

                // Set user data in p4est quadrant
                quadrant->p.user_data = (void*) newElement;
                // quadrant->p.user_data = (void*) &(elements[ELEMENT_COUNTER]);
                // quadrant->p.user_int = ELEMENT_COUNTER;
                // std::cout << "quadrant->p.user_data = " << quadrant->p.user_data << std::endl;

                // Update counter
                ELEMENT_COUNTER++;

                return;

            },
            this
        );

        nElements_ = ELEMENT_COUNTER;
        
        // Write initial mesh to file
        p4est_vtk_write_file(p4est_, NULL, "initial_mesh");

        // Set element neighbors
        p4est_iterate(
            p4est_,
            NULL,
            NULL,
            // [](p4est_iter_volume_info_t* info, void* user_data){
            //     p4est_t* p4est = info->p4est;
            //     p4est_quadrant_t* quadrant = info->quad;
            //     ElementMesh2D<FloatingDataType>& mesh = *((ElementMesh2D<FloatingDataType>*) p4est->user_pointer);

            //     // std::cout << "quadrant ID = " << quadrant->p.user_int;
            //     // std::cout << "quadrant->p.user_data = " << quadrant->p.user_data << std::endl;
            //     QuadElement2D<FloatingDataType>& element = *(static_cast<QuadElement2D<FloatingDataType>*>(quadrant->p.user_data));
            //     // QuadElement2D<FloatingDataType>& element = mesh.elements()[quadrant->p.user_int];

            //     // std::cout << "quadrant->p.user_int = " << quadrant->p.user_int << std::endl;
            //     // std::cout << "element = " << element.ID() << std::endl;
            //     // for (auto i = 0; i < 4; i++) std::cout << "  p" << i << " = [" << element.points()[i].x() << ", " << element.points()[i].y() << "]" << std::endl;

            //     return;
            // },
            NULL,
            [](p4est_iter_face_info_t* info, void* user_data) {
                p4est_t* p4est = info->p4est;
                sc_array_t* sides = &(info->sides);
                p4est_iter_face_side_t* side[2];
                ElementMesh2D<FloatingDataType>& mesh = *((ElementMesh2D<FloatingDataType>*) p4est->user_pointer);

                // std::cout << "sides->elem_count = " << sides->elem_count << std::endl;
                side[0] = p4est_iter_fside_array_index_int(sides, 0);
                side[1] = p4est_iter_fside_array_index_int(sides, 1);

                // std::cout << "side[0]->face = " << (int) side[0]->face << std::endl;
                // std::cout << "side[1]->face = " << (int) side[1]->face << std::endl;

                int p4estFaceIndex = side[0]->face;
                int internalFaceIndex = p4est::p4est_to_HydroForest_face_index(side[0]->face);
                int externalFaceIndex = p4est::p4est_to_HydroForest_face_index(side[1]->face);

                // std::cout << "p4estFaceIndex = " << p4estFaceIndex << std::endl;
                // std::cout << "faceIndex = " << faceIndex << std::endl;

                p4est_quadrant_t* quadrant;
                p4est_quadrant_t* neighborQuadrant; // THIS WILL NEED TO CHANGE FOR AMR

                // std::cout << "side[0]->is_hanging = " << (bool) side[0]->is_hanging << std::endl;
                // std::cout << "side[1]->is_hanging = " << (bool) side[1]->is_hanging << std::endl;

                if (side[0]->is_hanging == 0) {
                    // std::cout << "HERE" << std::endl;
                    quadrant = side[0]->is.full.quad;
                    neighborQuadrant = side[1]->is.full.quad;
                }

                // std::cout << "quadrant ID = " << quadrant->p.user_int;
                // std::cout << "neighbor ID = " << quadrant->p.user_int;

                QuadElement2D<FloatingDataType>& element = *((QuadElement2D<FloatingDataType>*) quadrant->p.user_data);
                QuadElement2D<FloatingDataType>& neighbor = *((QuadElement2D<FloatingDataType>*) neighborQuadrant->p.user_data);

                // std::cout << "element  = " << element.ID() << " @ l = " << internalFaceIndex << " : " << neighbor.ID() << std::endl;
                // std::cout << "neighbor = " << neighbor.ID() << " @ l = " << externalFaceIndex << " : " << element.ID() << std::endl << std::endl;
                element.neighbors()[internalFaceIndex] = &neighbor;
                neighbor.neighbors()[externalFaceIndex] = &element;

                // InternalBoundary<FloatingDataType> elementBoundary(neighbor);
                // InternalBoundary<FloatingDataType> neighborBoundary(element);
                // element.setBoundary(internalFaceIndex, elementBoundary);
                // neighbor.setBoundary(externalFaceIndex, neighborBoundary);


                // std::cout << "element = " << element.ID() << std::endl;
                // for (auto i = 0; i < 4; i++) std::cout << "  p" << i << " = [" << element.points()[i].x() << ", " << element.points()[i].y() << "]" << std::endl;

                // std::cout << "neighbor = " << neighbor.ID() << std::endl;
                // for (auto i = 0; i < 4; i++) std::cout << "  p" << i << " = [" << neighbor.points()[i].x() << ", " << neighbor.points()[i].y() << "]" << std::endl;

                return;

            },
            // NULL,
            NULL
        );

    }

    // ~ElementMesh2D() {
    //     // TODO
    //     // for (auto e = 0; e < elements_.size(); e++) {
    //     //     delete &elements_[e];
    //     // }
    // }

    FloatingDataType xLower() { return xLower_; }
    FloatingDataType xUpper() { return xUpper_; }
    FloatingDataType yLower() { return yLower_; }
    FloatingDataType yUpper() { return yUpper_; }
    int size() { return nElements_; }
    int order() { return order_; }
    // std::vector<QuadElement2D<FloatingDataType>>& elements() { return elements_; }

    void iterate(std::function<void(QuadElement2D<FloatingDataType>&)> elementCallbackFunction) {
        p4est_iterate(
            p4est_,
            NULL,
            &elementCallbackFunction,
            [](p4est_iter_volume_info_t* info, void* user_data) {
                // p4est_t* p4est = info->p4est;
                // ElementMesh2D<FloatingDataType>& mesh = *((ElementMesh2D<FloatingDataType>*) p4est->user_pointer);

                p4est_quadrant_t* quadrant = info->quad;
                QuadElement2D<FloatingDataType>& element = *(static_cast<QuadElement2D<FloatingDataType>*>(quadrant->p.user_data));
                std::function<void(QuadElement2D<FloatingDataType>&)>& elementCallbackFunction = *(static_cast<std::function<void(QuadElement2D<FloatingDataType>&)>*>(user_data));
                elementCallbackFunction(element);

                return;
            },
            NULL,
            NULL
        );
    }

    QuadElement2D<FloatingDataType>& element(int ID) {
        typedef struct element_context {
            int ID;
            QuadElement2D<FloatingDataType>* elementPtr;
        }
        element_context_t;

        element_context_t context;
        context.ID = ID;
        // QuadElement2D<FloatingDataType>* elementPtr;
        p4est_iterate(
            p4est_,
            NULL,
            &context,
            [](p4est_iter_volume_info_t* info, void* user_data) {
                p4est_quadrant_t* quadrant = info->quad;
                element_context_t* context = static_cast<element_context_t*>(user_data);
                QuadElement2D<FloatingDataType>& element = *(static_cast<QuadElement2D<FloatingDataType>*>(quadrant->p.user_data));
                if (element.ID() == context->ID) {
                    context->elementPtr = &element;
                }
                return;
            },
            NULL,
            NULL
        );
        return *context.elementPtr;
    }

    void setVariables(std::vector<std::string> variableNames) {
        iterate([&](QuadElement2D<FloatingDataType>& element){
            for (auto s = 0; s < variableNames.size(); s++) {
                element.vectorMap()[variableNames[s]] = Vector<FloatingDataType>(element.referenceGrid().size(), 0);
                element.vector(variableNames[s]).setName(variableNames[s]);
            }
        });
        // for (auto e = 0; e < nElements_; e++) {
        //     QuadElement2D<FloatingDataType>& element = elements_[e];
        //     for (auto s = 0; s < variableNames.size(); s++) {
        //         element.vectorMap()[variableNames[s]] = Vector<FloatingDataType>(element.referenceGrid().size(), 0);
        //         element.vector(variableNames[s]).setName(variableNames[s]);
        //     }
        // }
    }

    void setSolution(std::string variableName, std::function<FloatingDataType(FloatingDataType, FloatingDataType)> function) {
        iterate([&](QuadElement2D<FloatingDataType>& element){
            LobattoTensorProductGrid2D<FloatingDataType>& grid = element.referenceGrid();
            for (auto i = 0; i < grid.size(); i++) {
                SpaceVector2D<FloatingDataType> xi = grid(i);
                SpaceVector2D<FloatingDataType> x = element.mapReference2Physical(xi);
                element.vectorMap()[variableName][i] = function(x.x(), x.y());
            }
        });
        // for (auto e = 0; e < nElements_; e++) {
        //     QuadElement2D<FloatingDataType>& element = elements_[e];
        //     LobattoTensorProductGrid2D<FloatingDataType>& grid = element.referenceGrid();
        //     for (auto i = 0; i < grid.size(); i++) {
        //         SpaceVector2D<FloatingDataType> xi = grid(i);
        //         SpaceVector2D<FloatingDataType> x = element.mapReference2Physical(xi);
        //         element.vectorMap()[variableName][i] = function(x.x(), x.y());
        //     }
        // }

    }

    // void setMassMatrices() {
    //     for (auto& element : elements_) {
    //         element.matrixMap()["mass"] = DGMassMatrix2D<FloatingDataType>(element);
    //     }
    // }

    // void setDerivativeMatrices() {
    //     for (auto& element : elements_) {
    //         element.matrixMap()["derivative-x"] = DGDerivativeMatrix2D<FloatingDataType>(element, 0);
    //         element.matrixMap()["derivative-y"] = DGDerivativeMatrix2D<FloatingDataType>(element, 1);
    //     }
    // }

    void toVTK(std::string filename, std::vector<std::string> variableNames) {

        // Create header file
        XMLNode nodeVTKFile("VTKFile");
        nodeVTKFile.addAttribute("type", "PRectilinearGrid");
        nodeVTKFile.addAttribute("version", "0.1");
        nodeVTKFile.addAttribute("byte_order", "LittleEndian");

        XMLNode nodePRectilinearGrid("PRectilinearGrid");
        nodePRectilinearGrid.addAttribute("WholeExtent", "0 " + std::to_string(sqrt(nElements_)*order_) + " 0 " + std::to_string(sqrt(nElements_)*order_) + " 0 0");
        nodePRectilinearGrid.addAttribute("GhostLevel", "0");

        XMLNode nodePPointData("PPointData");
        std::vector<XMLNode> nodesPPointDataArray;
        std::string names = "";
        for (auto s = 0; s < variableNames.size(); s++) {
            names += variableNames[s] + " ";
            nodesPPointDataArray.push_back(XMLNode("PDataArray"));
            nodesPPointDataArray[s].addAttribute("type", "Float32");
            nodesPPointDataArray[s].addAttribute("Name", variableNames[s]);
            nodesPPointDataArray[s].addAttribute("NumberOfComponents", "1");
        }
        names += "elementID";
        nodePPointData.addAttribute("Scalars", names);

        XMLNode nodePCoordinates("PCoordinates");
        XMLNode nodeXCoordDataArray("PDataArray");
        nodeXCoordDataArray.addAttribute("type", "Float32");
        nodeXCoordDataArray.addAttribute("NumberOfComponents", "1");
        XMLNode nodeYCoordDataArray("PDataArray");
        nodeYCoordDataArray.addAttribute("type", "Float32");
        nodeYCoordDataArray.addAttribute("NumberOfComponents", "1");
        XMLNode nodeZCoordDataArray("PDataArray");
        nodeZCoordDataArray.addAttribute("type", "Float32");
        nodeZCoordDataArray.addAttribute("NumberOfComponents", "1");

        std::vector<XMLNode> nodesPiece;

        int e = 0;
        iterate([&](QuadElement2D<FloatingDataType>& element){
            LobattoTensorProductGrid2D<FloatingDataType>& grid = element.physicalGrid();

            RectilinearGridVTK gridVTK;
            gridVTK.buildMesh(grid, &grid.xGrid(), &grid.yGrid(), nullptr);
            for (auto s = 0; s < variableNames.size(); s++) {
                // element.vector(variableNames[s]).setName(variableNames[s]);
                gridVTK.addPointData(element.vector(variableNames[s]));
            }
            Vector<int> elementIDVector(element.size(), element.ID());
            elementIDVector.setName("elementID");
            gridVTK.addPointData(elementIDVector);
            std::string completeFilename = filename + "_" + std::to_string(e) + ".vtr";
            gridVTK.toVTK(completeFilename);

            nodesPiece.push_back(XMLNode("Piece"));
            nodesPiece[e].addAttribute("Extent", "0 " + std::to_string(order_) + " 0 " + std::to_string(order_) + " 0 0");
            nodesPiece[e].addAttribute("Source", completeFilename);

            e++;
        });
        // for (auto e = 0; e < nElements_; e++) {
        //     QuadElement2D<FloatingDataType>& element = elements_[e];
        //     LobattoTensorProductGrid2D<FloatingDataType>& grid = element.physicalGrid();

        //     RectilinearGridVTK gridVTK;
        //     gridVTK.buildMesh(grid, &grid.xGrid(), &grid.yGrid(), nullptr);
        //     for (auto s = 0; s < variableNames.size(); s++) {
        //         // element.vector(variableNames[s]).setName(variableNames[s]);
        //         gridVTK.addPointData(element.vector(variableNames[s]));
        //     }
        //     Vector<int> elementIDVector(element.size(), element.ID());
        //     elementIDVector.setName("elementID");
        //     gridVTK.addPointData(elementIDVector);
        //     std::string completeFilename = filename + "_" + std::to_string(e) + ".vtr";
        //     gridVTK.toVTK(completeFilename);

        //     nodesPiece.push_back(XMLNode("Piece"));
        //     nodesPiece[e].addAttribute("Extent", "0 " + std::to_string(order_) + " 0 " + std::to_string(order_) + " 0 0");
        //     nodesPiece[e].addAttribute("Source", completeFilename);

        // }

        nodeVTKFile.addChild(nodePRectilinearGrid);
            nodePRectilinearGrid.addChild(nodePPointData);
                for (auto& node : nodesPPointDataArray) nodePPointData.addChild(node);
            nodePRectilinearGrid.addChild(nodePCoordinates);
                nodePCoordinates.addChild(nodeXCoordDataArray);
                nodePCoordinates.addChild(nodeYCoordDataArray);
                nodePCoordinates.addChild(nodeZCoordDataArray);
            for (auto& node : nodesPiece) nodePRectilinearGrid.addChild(node);

        XMLTree tree(nodeVTKFile);
        tree.write(filename + ".pvtr");
    }

private:

    // static void p4est_element_init_(p4est_t* p4est, p4est_topidx_t which_tree, p4est_quadrant_t* quadrant) {
        
    //     // Get mesh and elements
    //     ElementMesh2D<FloatingDataType>& mesh = *((ElementMesh2D<FloatingDataType>*) p4est->user_pointer);
    //     std::vector<QuadElement2D<FloatingDataType>>& elements = mesh.elements();

    //     // Get quadrant points
    //     p4est_qcoord_t quad_length = P4EST_QUADRANT_LEN(quadrant->level);
    //     double xyz0[3];
    //     double xyz1[3];
    //     double xyz2[3];
    //     double xyz3[3];
    //     p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x, quadrant->y, xyz0);
    //     p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x + quad_length, quadrant->y, xyz1);
    //     p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x + quad_length, quadrant->y + quad_length, xyz2);
    //     p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x, quadrant->y + quad_length, xyz3);
        
    //     // Create vector of points
    //     std::vector<SpaceVector2D<FloatingDataType>> points = {
    //         SpaceVector2D<FloatingDataType>(xyz0[0], xyz0[1]),
    //         SpaceVector2D<FloatingDataType>(xyz1[0], xyz1[1]),
    //         SpaceVector2D<FloatingDataType>(xyz2[0], xyz2[1]),
    //         SpaceVector2D<FloatingDataType>(xyz3[0], xyz3[1])
    //     };

    //     // Create element and add to elements
    //     elements.push_back(QuadElement2D<FloatingDataType>(points, mesh.order(), mesh.order(), ELEMENT_COUNTER));

    //     // Set user data in p4est quadrant
    //     quadrant->p.user_data = &elements[ELEMENT_COUNTER];
    //     quadrant->p.user_int = ELEMENT_COUNTER;

    //     // Update counter
    //     ELEMENT_COUNTER++;

    // }

};

} // NAMESPACE : HydroForest

#endif // MESH_2D_HPP_