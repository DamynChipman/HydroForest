#include <iostream>
#include <HydroForestApp.hpp>
#include <UniformGrid1D.hpp>
#include <Options.hpp>
#include <XMLTree.hpp>

int main(int argc, char** argv) {

    HydroForest::HydroForestApp app(&argc, &argv);

    HydroForest::UniformGrid1D<double> xGrid(-1, 1, 10);
    HydroForest::UniformGrid1D<double> yGrid(2, 5, 10);
    HydroForest::UniformGrid1D<double> zGrid(0, 0, 0);
    HydroForest::RectilinearGridVTK vtk;
    vtk.buildMesh(xGrid, xGrid, yGrid, zGrid);
    vtk.toVTK("mesh.vtr");

    // HydroForest::XMLNode VTKNode("VTKFile");
    // VTKNode.addAttribute("type", "PolyData");
    // VTKNode.addAttribute("version", "0.1");
    // VTKNode.addAttribute("byte_order", "LittleEndian");

    // HydroForest::XMLNode PolyDataNode("PolyData");

    // HydroForest::XMLNode PieceNode("Piece");
    // PieceNode.addAttribute("NumberOfPoints", "8");
    // PieceNode.addAttribute("NumberOfVerts", "0");
    // PieceNode.addAttribute("NumberOfLines", "0");
    // PieceNode.addAttribute("NumberOfStrips", "0");
    // PieceNode.addAttribute("NumberOfPolys", "6");

    // HydroForest::XMLNode PointsNode("Points");
    
    // HydroForest::XMLNode PointsDataArrayNode("DataArray");
    // PointsDataArrayNode.addAttribute("type", "Float32");
    // PointsDataArrayNode.addAttribute("NumberOfComponents", "3");
    // PointsDataArrayNode.addAttribute("format", "ascii");
    // PointsDataArrayNode.data = "0 0 0 1 0 0 1 1 0 0 1 0 0 0 1 1 0 1 1 1 1 0 1 1";

    // HydroForest::XMLNode PointDataNode("PointData");
    // PointDataNode.addAttribute("Scalars", "my_scalars");

    // HydroForest::XMLNode PointDataDataArrayNode("DataArray");
    // PointDataDataArrayNode.addAttribute("type", "Float32");
    // PointDataDataArrayNode.addAttribute("Name", "my_scalars");
    // PointDataDataArrayNode.addAttribute("format", "ascii");
    // PointDataDataArrayNode.data = "0 1 2 3 4 5 6 7";

    // HydroForest::XMLNode CellDataNode("CellData");
    // CellDataNode.addAttribute("Scalars", "cell_scalars");
    // CellDataNode.addAttribute("Normals", "cell_normals");

    // HydroForest::XMLNode CellDataDataArrayNode1("DataArray");
    // CellDataDataArrayNode1.addAttribute("type", "Int32");
    // CellDataDataArrayNode1.addAttribute("Name", "cell_scalars");
    // CellDataDataArrayNode1.addAttribute("format", "ascii");
    // CellDataDataArrayNode1.data = "0 1 2 3 4 5";

    // HydroForest::XMLNode CellDataDataArrayNode2("DataArray");
    // CellDataDataArrayNode2.addAttribute("type", "Float32");
    // CellDataDataArrayNode2.addAttribute("Name", "cell_normals");
    // CellDataDataArrayNode2.addAttribute("NumberOfComponents", "3");
    // CellDataDataArrayNode2.addAttribute("format", "ascii");
    // CellDataDataArrayNode2.data = "0 0 -1 0 0 1 0 -1 0 0 1 0 -1 0 0 1 0 0";

    // HydroForest::XMLNode PolysNode("Polys");

    // HydroForest::XMLNode PolysDataArrayNode1("DataArray");
    // PolysDataArrayNode1.addAttribute("type", "Int32");
    // PolysDataArrayNode1.addAttribute("Name", "connectivity");
    // PolysDataArrayNode1.addAttribute("format", "ascii");
    // PolysDataArrayNode1.data = "0 1 2 3 4 5 6 7 0 1 5 4 2 3 7 6 0 4 7 3 1 2 6 5";

    // HydroForest::XMLNode PolysDataArrayNode2("DataArray");
    // PolysDataArrayNode2.addAttribute("type", "Int32");
    // PolysDataArrayNode2.addAttribute("Name", "offsets");
    // PolysDataArrayNode2.addAttribute("format", "ascii");
    // PolysDataArrayNode2.data = "4 8 12 16 20 24";

    // VTKNode.addChild(PolyDataNode);
    //     PolyDataNode.addChild(PieceNode);
    //         PieceNode.addChild(PointsNode);
    //             PointsNode.addChild(PointsDataArrayNode);
    //         PieceNode.addChild(PointDataNode);
    //             PointDataNode.addChild(PointDataDataArrayNode);
    //         PieceNode.addChild(CellDataNode);
    //             CellDataNode.addChild(CellDataDataArrayNode1);
    //             CellDataNode.addChild(CellDataDataArrayNode2);
    //         PieceNode.addChild(PolysNode);
    //             PolysNode.addChild(PolysDataArrayNode1);
    //             PolysNode.addChild(PolysDataArrayNode2);

    // HydroForest::XMLTree xmlTree(VTKNode);
    // xmlTree.write("test.vtp");

    return EXIT_SUCCESS;
}