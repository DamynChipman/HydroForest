#include "VTKBuilder.hpp"

namespace HydroForest {

RectilinearGridVTK::RectilinearGridVTK() : meshComplete_(false), coordsDataVector_(0), pointDataVector_(0), cellDataVector_(0), root_("VTKFile") {
    root_.addAttribute("type", "RectilinearGrid");
    root_.addAttribute("version", "0.1");
    root_.addAttribute("byte_order", "LittleEndian");
}

void RectilinearGridVTK::buildMesh(RectilinearGridNodeBase& mesh, DataArrayNodeBase& xCoords, DataArrayNodeBase& yCoords, DataArrayNodeBase& zCoords) {
    mesh_ = &mesh;
    coordsDataVector_.push_back(&xCoords);
    coordsDataVector_.push_back(&yCoords);
    coordsDataVector_.push_back(&zCoords);
    meshComplete_ = true;
    return;
}

void RectilinearGridVTK::toVTK(std::string filename) {

    if (meshComplete_) {

        XMLNode nodeRectilinearGrid("RectilinearGrid");
        nodeRectilinearGrid.addAttribute("WholeExtent", mesh_->getWholeExtent());

        XMLNode nodePiece("Piece");
        nodePiece.addAttribute("Extent", mesh_->getExtent());

        XMLNode nodeCoordinates("Coordinates");
        
        root_.addChild(nodeRectilinearGrid);
        nodeRectilinearGrid.addChild(nodePiece);
        nodePiece.addChild(nodeCoordinates);
        for (auto i = 0; i < coordsDataVector_.size(); i++) nodeCoordinates.addChild(coordsDataVector_[i]->toVTK());

        XMLTree tree(root_);
        tree.write(filename);

    }
    else {
        throw std::invalid_argument("[RectilinearGridVTK::toVTK] mesh and data not complete.");
    }

}

XMLNode* DataArrayNodeBase::toVTK() {
    XMLNode* nodeDataArray = new XMLNode("DataArray");
    nodeDataArray->addAttribute("type", getType());
    nodeDataArray->addAttribute("Name", getName());
    nodeDataArray->addAttribute("NumberOfComponents", getNumberOfComponents());
    nodeDataArray->addAttribute("format", getFormat());
    nodeDataArray->addAttribute("RangeMin", getRangeMin());
    nodeDataArray->addAttribute("RangeMax", getRangeMax());
    nodeDataArray->data = getData();
    return nodeDataArray;
}

} // NAMESPACE: HydroForest