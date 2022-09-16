#ifndef UNIFORM_GRID_1D_HPP_
#define UNIFORM_GRID_1D_HPP_

#include "Grid1DBase.hpp"
#include "VTKBuilder.hpp"

namespace HydroForest {

template<typename DataType>
class UniformGrid1D : public Grid1DBase<DataType>, public RectilinearGridNodeBase, public DataArrayNodeBase {

public:

    UniformGrid1D(DataType lowerBound, DataType upperBound, std::size_t nPoints) : lowerBound_(lowerBound), upperBound_(upperBound), nPoints_(nPoints), points_(nPoints) {
        points_.resize(nPoints_);
        gridSpacing_ = (upperBound_ - lowerBound_) / (nPoints_ - 1);
        for (auto i = 0; i < nPoints_; i++) {
            points_[i] = lowerBound_ + i * gridSpacing_;
        }
    }

    std::vector<DataType> getPoints() { return points_; }
    std::size_t getNPoints() { return points_.size(); }
    DataType getLowerBound() { return lowerBound_; }
    DataType getUpperBound() { return upperBound_; }

    std::string getWholeExtent() {
        return "0 " + std::to_string(nPoints_) + " 0 0 0 0";
    }

    std::string getExtent() {
        return "0 " + std::to_string(nPoints_) + " 0 0 0 0";
    }

    std::string getType() {
        return "Float32";
    }

    std::string getName() {
        return "UniformGrid1d";
    }

    std::string getNumberOfComponents() {
        return "1";
    }

    std::string getFormat() {
        return "ascii";
    }

    std::string getRangeMin() {
        return std::to_string(lowerBound_);
    }

    std::string getRangeMax() {
        return std::to_string(upperBound_);
    }

    std::string getData() {
        if (nPoints_ == 0) {
            return "0";
        }
        else {
            std::string pointsAsString = "";
            for (auto& p : points_) pointsAsString += std::to_string(p) + " ";
            return pointsAsString;
        }
    }

private:

    DataType lowerBound_;
    DataType upperBound_;
    DataType gridSpacing_;
    std::size_t nPoints_;
    std::vector<DataType> points_;

};

} // NAMESPACE: HydroForest

#endif // UNIFORM_GRID_1D_HPP_