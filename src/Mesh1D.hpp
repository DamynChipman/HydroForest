#ifndef MESH_1D_HPP
#define MESH_1D_HPP

#include <string>
#include <matplotlibcpp.h>
#include "Element1D.hpp"

namespace plt = matplotlibcpp;

namespace HydroForest {

template<typename FloatingDataType>
class Mesh1DBase {



};

template<typename FloatingDataType>
class ElementMesh1D : public Mesh1DBase<FloatingDataType> {

public:

    ElementMesh1D(FloatingDataType xLower, FloatingDataType xUpper, std::size_t nElements, std::size_t order) :
        xLower_(xLower),
        xUpper_(xUpper),
        nElements_(nElements),
        order_(order),
        elements_() {

        // Create elements
        FloatingDataType dx = (xUpper_ - xLower_) / nElements_;
        for (auto i = 0; i < nElements_; i++) {
            FloatingDataType xL = xLower_ + i*dx;
            FloatingDataType xR = xL + dx;
            LobattoGrid1D<FloatingDataType>* nodalGrid = new LobattoGrid1D<FloatingDataType>(order);
            
            LobattoGrid1D<FloatingDataType>* quadratureGrid;
            HydroForestApp& app = HydroForestApp::getInstance();
            if (std::get<std::string>(app.getOptions()["integration"]) == "exact") {
                quadratureGrid = new LobattoGrid1D<FloatingDataType>(order+1);
            }
            else {
                quadratureGrid = new LobattoGrid1D<FloatingDataType>(order);
            }

            elements_.push_back(Element1D<FloatingDataType>(xL, xR, *nodalGrid, *quadratureGrid));
        }

    }

    ~ElementMesh1D() {
        for (auto& e : elements_) {
            delete e.grid();
            delete e.quadratureGrid();
        }
    }

    FloatingDataType xLower() const { return xLower_; }
    FloatingDataType xUpper() const { return xUpper_; }
    std::size_t nElements() const { return nElements_; }
    std::size_t size() const { return nElements_; }
    std::size_t order() const { return order_; }
    Element1D<FloatingDataType>& operator[](std::size_t index) { return elements_[index]; }
    std::vector<Element1D<FloatingDataType>>& elements() { return elements_; }
    const std::vector<Element1D<FloatingDataType>>& elements() const { return elements_; }

    void setInitialCondition(std::function<FloatingDataType(FloatingDataType)> function, std::function<FloatingDataType(FloatingDataType)> flux) {
        // Iterate over elements in mesh
        for (auto e = 0; e < nElements_; e++) {
            // Iterate over points in element
            for (auto i = 0; i < elements_[e].size(); i++) {
                Grid1DBase<FloatingDataType>* elementGrid = elements_[e].grid();
                FloatingDataType xi = elementGrid->operator[](i);
                FloatingDataType x = elements_[e].transformLocal2Global(xi);
                elements_[e][i] = function(x);
                elements_[e].flux()[i] = flux(elements_[e][i]);
            }
        }
        return;
    }

    void plot(std::string format) {
        for (auto& e : elements_) {
            std::vector<FloatingDataType> x = e.grid()->getPoints().data();
            std::vector<FloatingDataType> y = e.solution().data();
            for (auto i = 0; i < e.size(); i++) {
                x[i] = e.transformLocal2Global(x[i]);
            }
            plt::plot(x, y, format);

        }
        for (auto& e : elements_) {
            plt::plot({e.xLower(), e.xUpper()}, {0, 0}, "|-k");
        }
    }

    friend std::ostream& operator<<(std::ostream& os, ElementMesh1D& mesh) {
        os << "--=== MESH ===--" << std::endl;
        os << "Global: [" << mesh.xLower() << ", " << mesh.xUpper() << "]" << std::endl;
        os << "# of Elements = " << mesh.nElements() << ", Order = " << mesh.order() << std::endl;
        for (auto i = 0; i < mesh.elements().size(); i++) {
            os << "===> " << i << " " << mesh.elements()[i];
        }
        return os;
    }

protected:

    FloatingDataType xLower_;
    FloatingDataType xUpper_;
    std::size_t nElements_;
    std::size_t order_;
    std::vector<Element1D<FloatingDataType>> elements_;

};

} // NAMESPACE : HydroForest

#endif // MESH_1D_HPP