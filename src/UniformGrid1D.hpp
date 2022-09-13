#ifndef UNIFORM_GRID_1D_HPP_
#define UNIFORM_GRID_1D_HPP_

#include "Grid1DBase.hpp"

namespace HydroForest {

template<typename DataType>
class UniformGrid1D : public Grid1DBase<DataType> {

public:

    UniformGrid1D(DataType lowerBound, DataType upperBound, std::size_t nPoints) : lowerBound_(lowerBound), upperBound_(upperBound), nPoints_(nPoints), points_(nPoints) {
        points_.resize(nPoints_);
        gridSpacing_ = (upperBound_ - lowerBound_) / (nPoints_ - 1);
        for (auto i = 0; i < nPoints_; i++) {
            points_[i] = lowerBound_ + i * gridSpacing_;
        }
    }

    std::vector<DataType> getPoints() {
        return points_;
    }

    DataType getLowerBound() { return lowerBound_; }
    DataType getUpperBound() { return upperBound_; }

private:

    DataType lowerBound_;
    DataType upperBound_;
    DataType gridSpacing_;
    std::size_t nPoints_;
    std::vector<DataType> points_;

};

} // NAMESPACE: HydroForest

#endif // UNIFORM_GRID_1D_HPP_