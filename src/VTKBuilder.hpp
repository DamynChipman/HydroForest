#ifndef VTK_BUILDER_HPP_
#define VTK_BUILDER_HPP_

#include <map>
#include <vector>

#include "XMLTree.hpp"

namespace HydroForest {

enum PointCellDataTypeVTK {
    SCALARS,
    VECTORS,
    NORMALS,
    TENSORS,
    TCOORDS
};

class DataArrayNodeBase {
    
public:

    virtual std::string getType() = 0;
    virtual std::string getName() = 0;
    virtual std::string getNumberOfComponents() = 0;
    virtual std::string getFormat() = 0;
    virtual std::string getRangeMin() = 0;
    virtual std::string getRangeMax() = 0;
    virtual std::string getData() = 0;

    XMLNode* toVTK();

};

struct EmptyDataArrayNode : public DataArrayNodeBase {
    EmptyDataArrayNode();
    std::string getType();
    std::string getName();
    std::string getNumberOfComponents();
    std::string getFormat();
    std::string getRangeMin();
    std::string getRangeMax();
    std::string getData();
};

class RectilinearGridNodeBase {

public:

    virtual std::string getWholeExtent() = 0;
    virtual std::string getExtent() = 0;

};

class RectilinearGridVTK {

public:

    RectilinearGridVTK();
    void buildMesh(RectilinearGridNodeBase& mesh, DataArrayNodeBase* xCoords, DataArrayNodeBase* yCoords, DataArrayNodeBase* zCoords);
    // void buildMesh(RectilinearGridNodeBase& mesh, DataArrayNodeBase& xCoords);
    void addPointData(DataArrayNodeBase& pointData);
    void addCellData(DataArrayNodeBase& cellData);
    void toVTK(std::string filename);

private:

    bool meshComplete_;

    XMLNode root_;

    RectilinearGridNodeBase* mesh_;
    EmptyDataArrayNode emptyDataArray_{};
    std::vector<DataArrayNodeBase*> coordsDataVector_;
    std::vector<DataArrayNodeBase*> pointDataVector_;
    std::vector<DataArrayNodeBase*> cellDataVector_;

};

class PointsVTK {

public:

    PointsVTK();
    PointsVTK(std::vector<double> xPoints);
    PointsVTK(std::vector<double> xPoints, std::vector<double> yPoints);
    PointsVTK(std::vector<double> xPoints, std::vector<double> yPoints, std::vector<double> zPoints);

    std::string getType();
    std::string getName();
    std::string getNumberOfComponents();
    std::string getFormat();
    int getNumberOfPoints();
    std::string getData();

    XMLNode* toVTK();

protected:

    std::vector<double> xPoints_;
    std::vector<double> yPoints_;
    std::vector<double> zPoints_;

};

class CellsVTK {

public:

    CellsVTK();
    XMLNode* toVTK();

protected:

    DataArrayNodeBase* connectivity_;
    DataArrayNodeBase* offsets_;
    DataArrayNodeBase* types_;

};

class PointDataVTK {

public:

    PointDataVTK();
    
    void addEntry(PointCellDataTypeVTK dataType, DataArrayNodeBase& dataEntry);

    XMLNode* toVTK();

protected:

    std::map<PointCellDataTypeVTK, std::vector<DataArrayNodeBase*>> entries_;

};

class UnstructuredGridVTK {

public:

    UnstructuredGridVTK(PointsVTK& pointsVTK, PointDataVTK& pointDataVTK);
    void toVTK(std::string filename);

protected:

    bool meshComplete_ = false;

    XMLNode root_;

    PointsVTK& pointsVTK_;
    PointDataVTK& pointDataVTK_;

};

} // NAMESPACE: HydroForest

#endif // VTK_BUILDER_HPP_