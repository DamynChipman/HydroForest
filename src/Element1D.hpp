#ifndef ELEMENT_1D_HPP_
#define ELEMENT_1D_HPP_

#include "Grid1D.hpp"
#include "Polynomial.hpp"

namespace HydroForest {

template<typename FloatingDataType>
class Element1D {

public:

    Element1D() {}
    Element1D(Grid1DBase<FloatingDataType>& grid, Grid1DBase<FloatingDataType>& quadratureGrid) :
        nodalPointgrid_(grid),
        quadratureGrid_(quadratureGrid),
        lagrangePoly_(grid.getPoints()),
        solutionVector_(grid.getNPoints())
            {}

    Grid1DBase<FloatingDataType>& grid() { return nodalPointgrid_; }
    Grid1DBase<FloatingDataType>& quadratureGrid() { return quadratureGrid_; }
    LagrangePolynomial& polynomial() { return lagrangePoly_; }
    Vector<FloatingDataType>& solution() { return solutionVector_; }
    FloatingDataType& operator[](std::size_t index) { return solutionVector_[index]; }

protected:

    Grid1DBase<FloatingDataType>& nodalPointgrid_;
    Grid1DBase<FloatingDataType>& quadratureGrid_;
    LagrangePolynomial lagrangePoly_;
    Vector<FloatingDataType> solutionVector_;

};

template<typename FloatingDataType>
std::vector<Element1D<FloatingDataType>> createElementGrid(UniformGrid1D<FloatingDataType>& elementGrid, Grid1DBase<FloatingDataType>& basisGrid, Grid1DBase<FloatingDataType>& quadratureGrid) {
    
    std::size_t nElements = elementGrid.getNPoints();
    std::vector<Element1D<FloatingDataType>> elements;
    for (auto e = 0; e < nElements; e++) {
        elements.push_back(Element1D<FloatingDataType>(basisGrid, quadratureGrid));
    }
    return elements;

}

} // NAMESPACE : HydroForest

#endif // ELEMENT_1D_HPP_