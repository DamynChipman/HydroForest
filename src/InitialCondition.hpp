#ifndef INITIAL_CONDITION_HPP_
#define INITIAL_CONDITION_HPP_

#include <functional>
#include "Mesh1D.hpp"

namespace HydroForest {

template<typename FloatingDataType>
class InitialConditionBase {};

template<typename FloatingDataType>
class FunctionInitialCondition : public InitialConditionBase<FloatingDataType> {

public:

    FunctionInitialCondition() {}
    FunctionInitialCondition(std::function<FloatingDataType(FloatingDataType x)> function) :
        function_(function)
            {}

    void setInitialCondtionOnMesh(ElementMesh1D<FloatingDataType>& mesh) {
        // Iterate over elements in mesh
        for (auto e = 0; e < mesh.size(); e++) {
            // Iterate over points in element
            for (auto i = 0; i < mesh[e].size(); i++) {
                Grid1DBase<FloatingDataType>* elementGrid = mesh[e].grid();
                FloatingDataType xi = elementGrid->operator[](i);
                FloatingDataType x = mesh[e].transformLocal2Global(xi);
                mesh[e][i] = function_(x);
            }
        }
        return;
    }

protected:

    std::function<FloatingDataType(FloatingDataType x)> function_;

};

} // NAMESPACE : HydroForest

#endif // INITIAL_CONDITION_HPP_