#ifndef ELEMENT_2D_HPP_
#define ELEMENT_2D_HPP_

#include <vector>
#include <string>
#include <map>
#include "SpaceVector.hpp"
#include "Grid2D.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

namespace HydroForest {

template<typename FloatingDataType>
class QuadElement2D {

protected:

    std::vector<SpaceVector2D<FloatingDataType>> points_;
    std::vector<SpaceVector2D<FloatingDataType>> normals_;
    LobattoTensorProductGrid2D<FloatingDataType> referenceGrid_;
    LobattoTensorProductGrid2D<FloatingDataType> physicalGrid_;
    LobattoTensorProductGrid2D<FloatingDataType> quadratureGrid_;
    std::map<std::string, Vector<FloatingDataType>> vecs_;
    std::map<std::string, Matrix<FloatingDataType>> mats_;

public:

    QuadElement2D(std::vector<SpaceVector2D<FloatingDataType>> points, int xOrder, int yOrder) :
        points_(points), referenceGrid_(xOrder, yOrder), physicalGrid_(xOrder, yOrder), quadratureGrid_(xOrder+1, yOrder+1) {

        // Compute element normals
        for (auto f = 0; f < 4; f++) {
            int p1 = f;
            int p2 = (f + 1) % 4;
            FloatingDataType x1 = points_[p1].x();
            FloatingDataType y1 = points_[p1].y();
            FloatingDataType x2 = points_[p2].x();
            FloatingDataType y2 = points_[p2].y();
            FloatingDataType dx = x2 - x1;
            FloatingDataType dy = y2 - y1;
            normals_.push_back(SpaceVector2D<FloatingDataType>(dy, -dx));
            normals_[f].normalize();
        }

        // Create physical grid via mapping
        Vector<FloatingDataType>& physicalGridpointsX = physicalGrid_.xGrid().getPoints();
        Vector<FloatingDataType>& physicalGridpointsY = physicalGrid_.yGrid().getPoints();
        for (auto i = 0; i < physicalGridpointsX.size(); i++) {
            SpaceVector2D<FloatingDataType> xi(physicalGridpointsX[i], 0);
            SpaceVector2D<FloatingDataType> x = mapReference2Physical(xi);
            physicalGridpointsX[i] = x.x();
        }
        for (auto j = 0; j < physicalGridpointsX.size(); j++) {
            SpaceVector2D<FloatingDataType> xi(0, physicalGridpointsY[j]);
            SpaceVector2D<FloatingDataType> x = mapReference2Physical(xi);
            physicalGridpointsY[j] = x.y();
        }

        // Compute transformation matrices
        Matrix<FloatingDataType> dh_dxi = referenceGrid_.xPoly().derivative(quadratureGrid_.xPoints());
        Matrix<FloatingDataType> h_eta = referenceGrid_.yPoly()(quadratureGrid_.yPoints());
        Matrix<FloatingDataType> h_xi = referenceGrid_.xPoly()(quadratureGrid_.xPoints());
        Matrix<FloatingDataType> dh_deta = referenceGrid_.yPoly().derivative(quadratureGrid_.yPoints());
        Matrix<FloatingDataType> dpsi_dxi = kroneckerProduct(dh_dxi, h_eta);
        Matrix<FloatingDataType> dpsi_deta = kroneckerProduct(h_xi, dh_deta);

        // Compute metric terms (Algorithm 12.1)
        vecs_["dx_dxi"] = Vector<FloatingDataType>(quadratureGrid_.size(), 0);
        vecs_["dx_deta"] = Vector<FloatingDataType>(quadratureGrid_.size(), 0);
        vecs_["dy_dxi"] = Vector<FloatingDataType>(quadratureGrid_.size(), 0);
        vecs_["dy_deta"] = Vector<FloatingDataType>(quadratureGrid_.size(), 0);
        vecs_["jacobian"] = Vector<FloatingDataType>(quadratureGrid_.size(), 0);
        vecs_["dxi_dx"] = Vector<FloatingDataType>(quadratureGrid_.size(), 0);
        vecs_["dxi_dy"] = Vector<FloatingDataType>(quadratureGrid_.size(), 0);
        vecs_["deta_dx"] = Vector<FloatingDataType>(quadratureGrid_.size(), 0);
        vecs_["deta_dy"] = Vector<FloatingDataType>(quadratureGrid_.size(), 0);
        for (auto k = 0; k < quadratureGrid_.size(); k++) {
            for (auto i = 0; i < referenceGrid_.size(); i++) {
                SpaceVector2D<FloatingDataType> xi = referenceGrid_(i);
                SpaceVector2D<FloatingDataType> x = mapReference2Physical(xi);
                vecs_["dx_dxi"][k] += dpsi_dxi(i,k)*x.x();
                vecs_["dx_deta"][k] += dpsi_deta(i,k)*x.x();
                vecs_["dy_dxi"][k] += dpsi_dxi(i,k)*x.y();
                vecs_["dy_deta"][k] += dpsi_deta(i,k)*x.y();
            }
        }
        for (auto k = 0; k < quadratureGrid_.size(); k++) {
            // vecs_["jacobian"][k] = vecs_["dx_dxi"][k]*vecs_["dy_deta"][k] - vecs_["dx_deta"][k]*vecs_["dy_dxi"][k];
            vecs_["jacobian"][k] = vecs_["dx_deta"][k]*vecs_["dy_dxi"][k] - vecs_["dx_dxi"][k]*vecs_["dy_deta"][k];
            vecs_["dxi_dx"][k] = vecs_["dy_deta"][k] / vecs_["jacobian"][k];
            vecs_["dxi_dy"][k] = -vecs_["dx_deta"][k] / vecs_["jacobian"][k];
            vecs_["deta_dx"][k] = -vecs_["dy_dxi"][k] / vecs_["jacobian"][k];
            vecs_["deta_dy"][k] = vecs_["dx_dxi"][k] / vecs_["jacobian"][k];
        }

        // std::cout << "jacobian = " << vecs_["jacobian"] << std::endl;
        // std::cout << "dx_dxi = " << vecs_["dx_dxi"] << std::endl;
        // std::cout << "dx_deta = " << vecs_["dx_deta"] << std::endl;
        // std::cout << "dy_dxi = " << vecs_["dy_dxi"] << std::endl;
        // std::cout << "dy_deta = " << vecs_["dy_deta"] << std::endl;
        // std::cout << "dxi_dx = " << vecs_["dxi_dx"] << std::endl;
        // std::cout << "dxi_dy = " << vecs_["dxi_dy"] << std::endl;
        // std::cout << "deta_dx = " << vecs_["deta_dx"] << std::endl;
        // std::cout << "deta_dy = " << vecs_["deta_dy"] << std::endl;

    }

    std::vector<SpaceVector2D<FloatingDataType>> points() { return points_; }
    std::vector<SpaceVector2D<FloatingDataType>> normals() { return normals_; }
    LobattoTensorProductGrid2D<FloatingDataType>& referenceGrid() { return referenceGrid_; }
    LobattoTensorProductGrid2D<FloatingDataType>& quadratureGrid() { return quadratureGrid_; }
    LobattoTensorProductGrid2D<FloatingDataType>& physicalGrid() { return physicalGrid_; }
    std::map<std::string, Vector<FloatingDataType>>& vectorMap() { return vecs_; }
    std::map<std::string, Matrix<FloatingDataType>>& matrixMap() { return mats_; }
    Vector<FloatingDataType>& vector(std::string name) { return vecs_[name]; }
    Matrix<FloatingDataType>& matrix(std::string name) { return mats_[name]; }
    std::size_t size() { return referenceGrid_.size(); }

    SpaceVector2D<FloatingDataType> mapReference2Physical(SpaceVector2D<FloatingDataType> xi) {
        Vector<FloatingDataType> m({
            0.25*(1 - xi.x())*(1 - xi.y()),
            0.25*(1 + xi.x())*(1 - xi.y()),
            0.25*(1 + xi.x())*(1 + xi.y()),
            0.25*(1 - xi.x())*(1 + xi.y())
        });
        Vector<FloatingDataType> x({
            points_[0].x(),
            points_[1].x(),
            points_[2].x(),
            points_[3].x()
        });
        Vector<FloatingDataType> y({
            points_[0].y(),
            points_[1].y(),
            points_[2].y(),
            points_[3].y()
        });
        return {m*x, m*y};
    }

    SpaceVector2D<FloatingDataType> mapPhysical2Reference(SpaceVector2D<FloatingDataType> x) {
        return x;
    }

};

} // NAMESPACE : HydroForest

#endif // ELEMENT_2D_HPP_