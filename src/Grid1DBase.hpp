#ifndef GRID_1D_BASE_HPP_
#define GRID_1D_BASE_HPP_
#pragma once

#include <vector>

namespace HydroForest {

template<typename DataType>
class Grid1DBase {

public:

    virtual std::vector<DataType> getPoints() = 0;
    virtual DataType getLowerBound() = 0;
    virtual DataType getUpperBound() = 0;

};

} // NAMESPACE: HydroForest

#endif // GRID_1D_BASE_HPP_