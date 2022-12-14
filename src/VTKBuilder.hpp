#ifndef VTK_BUILDER_HPP_
#define VTK_BUILDER_HPP_

#include "XMLTree.hpp"

namespace HydroForest {

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

class RectilinearGridNodeBase {

public:

    virtual std::string getWholeExtent() = 0;
    virtual std::string getExtent() = 0;

};

class RectilinearGridVTK {

public:

    RectilinearGridVTK();
    void buildMesh(RectilinearGridNodeBase& mesh, DataArrayNodeBase& xCoords, DataArrayNodeBase& yCoords, DataArrayNodeBase& zCoords);
    // void buildMesh(RectilinearGridNodeBase& mesh, DataArrayNodeBase& xCoords);
    void addPointData(DataArrayNodeBase& pointData);
    void addCellData(DataArrayNodeBase& cellData);
    void toVTK(std::string filename);

private:

    bool meshComplete_;

    XMLNode root_;

    RectilinearGridNodeBase* mesh_;
    std::vector<DataArrayNodeBase*> coordsDataVector_;
    std::vector<DataArrayNodeBase*> pointDataVector_;
    std::vector<DataArrayNodeBase*> cellDataVector_;

};

} // NAMESPACE: HydroForest

#endif // VTK_BUILDER_HPP_