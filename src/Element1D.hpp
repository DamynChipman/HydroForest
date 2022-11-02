#ifndef ELEMENT_1D_HPP_
#define ELEMENT_1D_HPP_

#include <string>
#include <map>
#include "Grid1D.hpp"
#include "Polynomial.hpp"

namespace HydroForest {

template<typename FloatingDataType>
class Element1D {

public:

    Element1D() :
        lagrangePoly_({})
            {}

    Element1D(FloatingDataType xLower, FloatingDataType xUpper, Grid1DBase<FloatingDataType>& grid, Grid1DBase<FloatingDataType>& quadratureGrid) :
        xLower_(xLower),
        xUpper_(xUpper),
        dx_(xUpper - xLower),
        nodalPointgrid_(&grid),
        quadratureGrid_(&quadratureGrid),
        lagrangePoly_(grid.getPoints()),
        solutionVector_(grid.getNPoints()),
        fluxVector_(grid.getNPoints())
            {}

    // Element1D(Element1D& other) : 
    //     xLower_(other.xLower()),
    //     xUpper_(other.xUpper()),
    //     dx_(other.dx()),
    //     nodalPointgrid_(other.grid()),
    //     quadratureGrid_(other.quadratureGrid()),
    //     lagrangePoly_(other.polynomial()),
    //     solutionVector_(other.solution()),
    //     fluxVector_(other.flux())
    //         {}

    // Element1D(const Element1D& other) : 
    //     xLower_(other.xLower()),
    //     xUpper_(other.xUpper()),
    //     dx_(other.dx()),
    //     nodalPointgrid_(&other.grid()),
    //     quadratureGrid_(&other.quadratureGrid()),
    //     lagrangePoly_(other.polynomial()),
    //     solutionVector_(other.solution()),
    //     fluxVector_(other.flux())
    //         {}

    Element1D<FloatingDataType>& operator=(Element1D<FloatingDataType>& other) {
        if (&other != this) {
            xLower_ = other.xLower();
            xUpper_ = other.xUpper();
            dx_ = other.dx();
            nodalPointgrid_ = other.grid();
            quadratureGrid_ = other.quadratureGrid();
            lagrangePoly_ = other.polynomial();
            solutionVector_ = other.solution();
            fluxVector_ = other.flux();
            return *this;
        }
        return *this;
    }

    FloatingDataType xLower() const { return xLower_; }
    FloatingDataType xUpper() const { return xUpper_; }
    FloatingDataType dx() const { return dx_; }
    std::size_t size() const { return nodalPointgrid_->getNPoints(); }
    const Grid1DBase<FloatingDataType>* grid() const { return nodalPointgrid_; }
    Grid1DBase<FloatingDataType>* grid() { return nodalPointgrid_; }
    const Grid1DBase<FloatingDataType>* quadratureGrid() const { return quadratureGrid_; }
    Grid1DBase<FloatingDataType>* quadratureGrid() { return quadratureGrid_; }
    const LagrangePolynomial& polynomial() const { return lagrangePoly_; }
    LagrangePolynomial& polynomial() { return lagrangePoly_; }
    const Vector<FloatingDataType>& solution() const { return solutionVector_; }
    Vector<FloatingDataType>& solution() { return solutionVector_; }
    const Vector<FloatingDataType>& flux() const { return fluxVector_; }
    Vector<FloatingDataType>& flux() { return fluxVector_; }
    FloatingDataType& operator[](std::size_t index) { return solutionVector_[index]; }
    std::map<std::string, Vector<FloatingDataType>>& vectors() { return vectors_; }
    std::map<std::string, Matrix<FloatingDataType>>& matrices() { return matrices_; }

    FloatingDataType transformLocal2Global(FloatingDataType xi) {
        FloatingDataType xMid = (xLower_ + xUpper_) / 2.0;
        return (xMid + (0.5)*dx_*xi);
    }

    FloatingDataType transformGlobal2Local(FloatingDataType x) {
        return (2.0*x - (xLower_ + xUpper_)) / (xUpper_ - xLower_);
    }

    friend std::ostream& operator<<(std::ostream& os, Element1D<FloatingDataType>& element) {
        os << "Element:" << std::endl;
        os << "[" << element.xLower() << ", " << element.xUpper() << "] " << std::endl;
        os << "dx = " << element.dx() << std::endl;
        os << "Nodal Point Grid: " << std::endl;
        os << *element.grid();
        os << "Quadrature Grid: " << std::endl;
        os << *element.quadratureGrid();
        os << "Solution Vector: " << std::endl;
        for (auto i = 0; i < element.solution().size(); i++) {
            os << std::setprecision(4) << std::setw(12) << element[i];
            if (i % 10 == 9) os << std::endl;
        }
        os << std::endl;
        os << "Flux Vector: " << std::endl;
        for (auto i = 0; i < element.flux().size(); i++) {
            os << std::setprecision(4) << std::setw(12) << element.flux()[i];
            if (i % 10 == 9) os << std::endl;
        }
        os << std::endl;
        return os;
    }

protected:

    FloatingDataType xLower_;
    FloatingDataType xUpper_;
    FloatingDataType dx_;
    Grid1DBase<FloatingDataType>* nodalPointgrid_;
    Grid1DBase<FloatingDataType>* quadratureGrid_;
    LagrangePolynomial lagrangePoly_;
    Vector<FloatingDataType> solutionVector_;
    Vector<FloatingDataType> fluxVector_;
    std::map<std::string, Vector<FloatingDataType>> vectors_;
    std::map<std::string, Matrix<FloatingDataType>> matrices_;

};

// template<typename FloatingDataType>
// struct Element1DFactory {

//     Element1DFactory() {}

//     static Element1D<FloatingDataType>& newElement1D(FloatingDataType xLower, FloatingDataType xUpper, std::size_t order) {
        
//         // Create element from app options
//         HydroForestApp& app = HydroForestApp::getInstance();
//         Options& opts = app.getOptions();
//         std::string nodalGridTypeFlag = opts["nodal-grid-type"];
//         std::string quadratureGridTypeFlag = opts["quadrature-grid-type"];

//         Grid1DBase<FloatingDataType>* nodalGrid_p;
//         if (nodalGridTypeFlag == "Lobatto") {
//             nodalGrid_p = new LobattoGrid1D<FloatingDataType>(order);
//         }
//         else {
//             throw std::invalid_argument("[HydroForest::Element1DFactory::newElement1D] NOT IMPLEMENTED!");
//         }

//         Grid1DBase<FloatingDataType>* quadratureGrid_p;
//         if (quadratureGridTypeFlag == "Lobatto") {
//             quadratureGrid_p = new LobattoGrid1D<FloatingDataType>(order);
//         }
//         else {
//             throw std::invalid_argument("[HydroForest::Element1DFactory::newElement1D] NOT IMPLEMENTED!");
//         }

//         Element1D<FloatingDataType> element(xLower, xUpper, *nodalGrid_p, *quadratureGrid_p);
//         return element;

//     }

// }

// template<typename FloatingDataType>
// std::vector<Element1D<FloatingDataType>> createElementGrid(UniformGrid1D<FloatingDataType>& elementGrid, Grid1DBase<FloatingDataType>& basisGrid, Grid1DBase<FloatingDataType>& quadratureGrid) {
    
//     std::size_t nElements = elementGrid.getNPoints();
//     std::vector<Element1D<FloatingDataType>> elements;
//     for (auto e = 0; e < nElements; e++) {
//         elements.push_back(Element1D<FloatingDataType>(basisGrid, quadratureGrid));
//     }
//     return elements;

// }

} // NAMESPACE : HydroForest

#endif // ELEMENT_1D_HPP_