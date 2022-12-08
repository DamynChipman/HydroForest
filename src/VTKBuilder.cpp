#include "VTKBuilder.hpp"

namespace HydroForest {

EmptyDataArrayNode::EmptyDataArrayNode() {}
std::string EmptyDataArrayNode::getType() { return "Float32"; }
std::string EmptyDataArrayNode::getName() { return "Empty"; }
std::string EmptyDataArrayNode::getNumberOfComponents() { return "1"; }
std::string EmptyDataArrayNode::getFormat() { return "ascii"; }
std::string EmptyDataArrayNode::getRangeMin() { return "0.0"; }
std::string EmptyDataArrayNode::getRangeMax() { return "0.0"; }
std::string EmptyDataArrayNode::getData() { return "0.0"; }

RectilinearGridVTK::RectilinearGridVTK() : meshComplete_(false), coordsDataVector_(0), pointDataVector_(0), cellDataVector_(0), root_("VTKFile") {
    root_.addAttribute("type", "RectilinearGrid");
    root_.addAttribute("version", "0.1");
    root_.addAttribute("byte_order", "LittleEndian");
}

void RectilinearGridVTK::buildMesh(RectilinearGridNodeBase& mesh, DataArrayNodeBase* xCoords, DataArrayNodeBase* yCoords, DataArrayNodeBase* zCoords) {
    mesh_ = &mesh;
    
    if (xCoords != nullptr) {
        coordsDataVector_.push_back(xCoords);
    }
    else {
        coordsDataVector_.push_back(&emptyDataArray_);
    }

    if (yCoords != nullptr) {
        coordsDataVector_.push_back(yCoords);
    }
    else {
        coordsDataVector_.push_back(&emptyDataArray_);
    }

    if (zCoords != nullptr) {
        coordsDataVector_.push_back(zCoords);
    }
    else {
        coordsDataVector_.push_back(&emptyDataArray_);
    }

    meshComplete_ = true;
    return;
}

void RectilinearGridVTK::addPointData(DataArrayNodeBase& pointData) {
    pointDataVector_.push_back(&pointData);
}

void RectilinearGridVTK::addCellData(DataArrayNodeBase& cellData) {
    cellDataVector_.push_back(&cellData);
}

void RectilinearGridVTK::toVTK(std::string filename) {

    if (meshComplete_) {

        XMLNode nodeRectilinearGrid("RectilinearGrid");
        nodeRectilinearGrid.addAttribute("WholeExtent", mesh_->getWholeExtent());

        XMLNode nodePiece("Piece");
        nodePiece.addAttribute("Extent", mesh_->getExtent());

        XMLNode nodeCellData("CellData");
        if (!cellDataVector_.empty()) {
            std::string names = "";
            for (auto& cellData : cellDataVector_) names += cellData->getName() + " ";
            nodeCellData.addAttribute("Scalars", names);
        }
        // for (auto cellData : cellDataVector_) {
        //     nodeCellData.addAttribute("Scalars", cellData->getName());
        // }

        XMLNode nodePointData("PointData");
        if (!pointDataVector_.empty()) {
            std::string names = "";
            for (auto& pointData : pointDataVector_) names += pointData->getName() + " ";
            nodePointData.addAttribute("Scalars", names);
        }
        // for (auto pointData : pointDataVector_) {
        //     nodePointData.addAttribute("Scalars", pointData->getName());
        // }

        XMLNode nodeCoordinates("Coordinates");
        
        root_.addChild(nodeRectilinearGrid);
        nodeRectilinearGrid.addChild(nodePiece);
        if (!cellDataVector_.empty()) nodePiece.addChild(nodeCellData);
        if (!pointDataVector_.empty()) nodePiece.addChild(nodePointData);
        nodePiece.addChild(nodeCoordinates);
        for (auto cellData : cellDataVector_) nodeCellData.addChild(cellData->toVTK());
        for (auto pointData : pointDataVector_) nodePointData.addChild(pointData->toVTK());
        for (auto coordData : coordsDataVector_) nodeCoordinates.addChild(coordData->toVTK());
        // for (auto i = 0; i < coordsDataVector_.size(); i++) nodeCoordinates.addChild(coordsDataVector_[i]->toVTK());

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

PointsVTK::PointsVTK() :
    xPoints_{}, yPoints_{}, zPoints_{}
        {}

PointsVTK::PointsVTK(std::vector<double> xPoints) :
    xPoints_(xPoints), yPoints_{}, zPoints_{}
        {}

PointsVTK::PointsVTK(std::vector<double> xPoints, std::vector<double> yPoints) :
    xPoints_(xPoints), yPoints_(yPoints), zPoints_{}
        {}

PointsVTK::PointsVTK(std::vector<double> xPoints, std::vector<double> yPoints, std::vector<double> zPoints) :
    xPoints_(xPoints), yPoints_(yPoints), zPoints_(zPoints)
        {}

std::string PointsVTK::getType() { return "Float32"; }
std::string PointsVTK::getName() { return "Points"; }
std::string PointsVTK::getNumberOfComponents() { return "3"; }
std::string PointsVTK::getFormat() { return "ascii"; }
int PointsVTK::getNumberOfPoints() {
    if (xPoints_.empty()) {
        xPoints_.push_back(0);
    }
    if (yPoints_.empty()) {
        yPoints_.push_back(0);
    }
    if (zPoints_.empty()) {
        zPoints_.push_back(0);
    }
    return xPoints_.size()*yPoints_.size()*zPoints_.size();
}
std::string PointsVTK::getData() {
    std::string str = "";
                
    if (xPoints_.empty()) {
        xPoints_.push_back(0);
    }
    if (yPoints_.empty()) {
        yPoints_.push_back(0);
    }
    if (zPoints_.empty()) {
        zPoints_.push_back(0);
    }

    for (auto i = 0; i < xPoints_.size(); i++) {
        for (auto j = 0; j < yPoints_.size(); j++) {
            for (auto k = 0; k < zPoints_.size(); k++) {
                str += std::to_string(xPoints_[i]) + " " + std::to_string(yPoints_[j]) + " " + std::to_string(zPoints_[k]) + "\n";
            }
        }
    }

    str.pop_back();
    return str;
}

XMLNode* PointsVTK::toVTK() {
    XMLNode* nodePoints = new XMLNode("Points");
    XMLNode* nodeDataArray = new XMLNode("DataArray");
    nodeDataArray->addAttribute("type", getType());
    nodeDataArray->addAttribute("Name", getName());
    nodeDataArray->addAttribute("NumberOfComponents", getNumberOfComponents());
    nodeDataArray->addAttribute("format", getFormat());
    nodeDataArray->data = getData();
    nodePoints->addChild(nodeDataArray);
    return nodePoints;
}

CellsVTK::CellsVTK() {}

XMLNode* CellsVTK::toVTK() {
    XMLNode* nodeCells = new XMLNode("Cells");
    
    XMLNode* nodeDataArrayConnectivity = new XMLNode("DataArray");
    nodeDataArrayConnectivity->addAttribute("type", "Int32");
    nodeDataArrayConnectivity->addAttribute("Name", "connectivity");
    nodeDataArrayConnectivity->addAttribute("NumberOfComponents", "1");
    nodeDataArrayConnectivity->addAttribute("format", "ascii");

    XMLNode* nodeDataArrayOffsets = new XMLNode("DataArray");
    nodeDataArrayOffsets->addAttribute("type", "Int32");
    nodeDataArrayOffsets->addAttribute("Name", "offsets");
    nodeDataArrayOffsets->addAttribute("NumberOfComponents", "1");
    nodeDataArrayOffsets->addAttribute("format", "ascii");

    XMLNode* nodeDataArrayTypes = new XMLNode("DataArray");
    nodeDataArrayTypes->addAttribute("type", "Int32");
    nodeDataArrayTypes->addAttribute("Name", "types");
    nodeDataArrayTypes->addAttribute("NumberOfComponents", "1");
    nodeDataArrayTypes->addAttribute("format", "ascii");

    nodeCells->addChild(nodeDataArrayConnectivity);
    nodeCells->addChild(nodeDataArrayOffsets);
    nodeCells->addChild(nodeDataArrayTypes);

    return nodeCells;
}

PointDataVTK::PointDataVTK() {}

void PointDataVTK::addEntry(PointCellDataTypeVTK dataType, DataArrayNodeBase& dataEntry) { entries_[dataType].push_back(&dataEntry); }

XMLNode* PointDataVTK::toVTK() {
    XMLNode* nodePointData = new XMLNode("PointData");
    std::string attributes;

    attributes = "";
    for (auto& dataArray : entries_[PointCellDataTypeVTK::SCALARS]) { attributes += dataArray->getName() + " "; }
    if (attributes != "") { nodePointData->addAttribute("Scalars", attributes); }

    attributes = "";
    for (auto& dataArray : entries_[PointCellDataTypeVTK::VECTORS]) { attributes += dataArray->getName() + " "; }
    if (attributes != "") { nodePointData->addAttribute("Vectors", attributes); }

    attributes = "";
    for (auto& dataArray : entries_[PointCellDataTypeVTK::NORMALS]) { attributes += dataArray->getName() + " "; }
    if (attributes != "") { nodePointData->addAttribute("Normals", attributes); }

    attributes = "";
    for (auto& dataArray : entries_[PointCellDataTypeVTK::TENSORS]) { attributes += dataArray->getName() + " "; }
    if (attributes != "") { nodePointData->addAttribute("Tensors", attributes); }

    attributes = "";
    for (auto& dataArray : entries_[PointCellDataTypeVTK::TCOORDS]) { attributes += dataArray->getName() + " "; }
    if (attributes != "") { nodePointData->addAttribute("TCoords", attributes); }

    for (auto& [type, vec] : entries_) {
        for (auto& dataArray : vec) {
            nodePointData->addChild(dataArray->toVTK());
        }
    }

    return nodePointData;
}

UnstructuredGridVTK::UnstructuredGridVTK(PointsVTK& pointsVTK, PointDataVTK& pointDataVTK) :
    root_("VTKFile"), pointsVTK_(pointsVTK), pointDataVTK_(pointDataVTK), meshComplete_(true)
        {
            root_.addAttribute("type", "UnstructuredGrid");
            root_.addAttribute("version", "0.1");
            root_.addAttribute("byte_order", "LittleEndian");
        }

void UnstructuredGridVTK::toVTK(std::string filename) {

    if (meshComplete_) {
        XMLNode nodeUnstructuredGrid("UnstructuredGrid");

        XMLNode nodePiece("Piece");
        nodePiece.addAttribute("NumberOfPoints", std::to_string(pointsVTK_.getNumberOfPoints()));
        nodePiece.addAttribute("NumberOfCells", "0");

        CellsVTK cellsVTK;

        root_.addChild(nodeUnstructuredGrid);
            nodeUnstructuredGrid.addChild(nodePiece);
                nodePiece.addChild(pointsVTK_.toVTK());
                nodePiece.addChild(cellsVTK.toVTK());
                nodePiece.addChild(pointDataVTK_.toVTK());

        XMLTree tree(root_);
        tree.write(filename);

    }

}
        

} // NAMESPACE: HydroForest