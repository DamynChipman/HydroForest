#ifndef GRID_2D_HPP_
#define GRID_2D_HPP_

#include <iterator>
#include <utility>
#include <cstddef>

#include "Grid1D.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "SpaceVector.hpp"
#include "Polynomial.hpp"
#include "VTKBuilder.hpp"

namespace HydroForest {

template<typename FloatingDataType>
class LobattoTensorProductGrid2D : public RectilinearGridNodeBase {

public:

    // LobattoTensorProductGrid2D() {}
    LobattoTensorProductGrid2D(std::size_t xOrder, std::size_t yOrder) :
        xGrid_(xOrder), yGrid_(yOrder), xPoly_(xGrid_.getPoints()), yPoly_(yGrid_.getPoints()), faceIndexSets_(4) {

        //
        faceIndexSets_[0] = Vector<int>(xSize()); // South
        faceIndexSets_[1] = Vector<int>(ySize()); // East
        faceIndexSets_[2] = Vector<int>(xSize()); // North
        faceIndexSets_[3] = Vector<int>(ySize()); // West

        for (auto i = 0; i < xSize(); i++) faceIndexSets_[0][i] = i;
        for (auto j = 0; j < ySize(); j++) faceIndexSets_[1][j] = (j + 1)*(xSize()) - 1;
        for (auto i = 0; i < xSize(); i++) faceIndexSets_[2][i] = (ySize() - 1)*ySize() + i;
        for (auto j = 0; j < ySize(); j++) faceIndexSets_[3][j] = j*xSize();

    }

    // struct Iterator
    // {
    //     using IteratorCategory = std::input_iterator_tag;
    //     using DifferenceType = std::ptrdiff_t;
    //     using ValueType = FloatingDataType;
    //     using Pointer = FloatingDataType*;
    //     using Reference = FloatingDataType&;

    //     Iterator(Pointer ptr) :
    //         ptr_(ptr)
    //             {}

    //     Reference operator*() const { return *ptr_; }
    //     Pointer operator->() { return ptr_; }
    //     Iterator& operator++() { ptr_++; return *this; }


    // private:

    //     Pointer ptr_;

    // };

    SpaceVector2D<FloatingDataType> operator()(std::size_t i, std::size_t j) {
        return {xGrid_[i], yGrid_[j]};
    }

    SpaceVector2D<FloatingDataType> operator()(std::size_t k) {
        return {xGrid_[ID(k).first], yGrid_[ID(k).second]};
    }

    LagrangePolynomial& xPoly() { return xPoly_; }
    LobattoGrid1D<FloatingDataType>& xGrid() { return xGrid_; }
    FloatingDataType xLower() { return xGrid_.getLowerBound(); }
    FloatingDataType xUpper() { return xGrid_.getUpperBound(); }
    Vector<FloatingDataType> xPoints() { return xGrid_.getPoints(); }
    Vector<FloatingDataType> xWeights() { return xGrid_.getWeights(); }
    std::size_t xOrder() { return xGrid_.order(); }
    std::size_t xSize() { return xGrid_.getNPoints(); }

    LagrangePolynomial& yPoly() { return yPoly_; }
    LobattoGrid1D<FloatingDataType>& yGrid() { return yGrid_; }
    FloatingDataType yLower() { return yGrid_.getLowerBound(); }
    FloatingDataType yUpper() { return yGrid_.getUpperBound(); }
    Vector<FloatingDataType> yPoints() { return yGrid_.getPoints(); }
    Vector<FloatingDataType> yWeights() { return yGrid_.getWeights(); }
    std::size_t yOrder() { return yGrid_.order(); }
    std::size_t ySize() { return yGrid_.getNPoints(); }

    std::size_t size() { return xGrid_.getNPoints() * yGrid_.getNPoints(); }
    std::vector<Vector<int>>& faceIndexSets() { return faceIndexSets_; }

    Matrix<FloatingDataType> basisMatrix(LobattoTensorProductGrid2D<FloatingDataType>& quadratureGrid) {
        Matrix<FloatingDataType> psi_x = xPoly_(quadratureGrid.xPoints());
        Matrix<FloatingDataType> psi_y = yPoly_(quadratureGrid.yPoints());
        return kroneckerProduct(psi_x, psi_y);
    }

    Matrix<FloatingDataType> basisDerivativeMatrix(LobattoTensorProductGrid2D<FloatingDataType>& quadratureGrid) {
        Matrix<FloatingDataType> psi_x = xPoly_.derivative(quadratureGrid.xPoints());
        Matrix<FloatingDataType> psi_y = yPoly_.derivative(quadratureGrid.yPoints());
        return kroneckerProduct(psi_x, psi_y);
    }

    Matrix<FloatingDataType> basisWeights() {
        return outerProduct(xWeights(), yWeights());
    }

    std::string getWholeExtent() {
        return "0 " + std::to_string(xGrid_.order()) + " 0 " + std::to_string(yGrid_.order()) + " 0 0";
    }

    std::string getExtent() {
        return "0 " + std::to_string(xGrid_.order()) + " 0 " + std::to_string(yGrid_.order()) + " 0 0";
    }

    int ID(int i, int j) { return j + i*ySize(); }
    std::pair<int, int> ID(int k) { return{k%xSize(), k/(int)ySize()}; }

protected:

    LobattoGrid1D<FloatingDataType> xGrid_;
    LobattoGrid1D<FloatingDataType> yGrid_;
    LagrangePolynomial xPoly_;
    LagrangePolynomial yPoly_;
    std::vector<Vector<int>> faceIndexSets_;

};

} // NAMESPACE : HydroForest

#endif // GRID_2D_HPP_