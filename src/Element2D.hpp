#ifndef ELEMENT_2D_HPP_
#define ELEMENT_2D_HPP_

#include <cmath>
#include <vector>
#include <string>
#include <map>
#include "SpaceVector.hpp"
// #include "Boundary2D.hpp"
#include "Grid2D.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

namespace HydroForest {

int neighborFaceIndex(int faceIndex) {
    switch (faceIndex) {
        case 0 : return 2;
        case 1 : return 3;
        case 2 : return 0;
        case 3 : return 1;
    }
}

// Forward declarations
template<typename FloatingDataType> class QuadElement2D;

template <typename FloatingDataType>
class BoundaryBase {
public:
    virtual QuadElement2D<FloatingDataType>& getElementAtBoundary() = 0;
};

template<typename FloatingDataType>
class InternalBoundary : public BoundaryBase<FloatingDataType> {

protected:

    QuadElement2D<FloatingDataType>& neighbor_;

public:

    InternalBoundary(QuadElement2D<FloatingDataType>& neighbor) :
        neighbor_(neighbor)
            {}

    QuadElement2D<FloatingDataType>& getElementAtBoundary() { return neighbor_; }

};

template<typename FloatingDataType>
class PeriodicBoundary : public BoundaryBase<FloatingDataType> {

protected:

    QuadElement2D<FloatingDataType>& neighbor_;

public:

    PeriodicBoundary(QuadElement2D<FloatingDataType>& neighbor) :
        neighbor_(neighbor)
            {}

    QuadElement2D<FloatingDataType>& getElementAtBoundary() { return neighbor_; }

};

template<typename FloatingDataType>
class QuadElement2D {

protected:

    int ID_;
    std::vector<SpaceVector2D<FloatingDataType>> points_;
    std::vector<SpaceVector2D<FloatingDataType>> normals_;
    std::vector<BoundaryBase<FloatingDataType>*> boundaries_;
    std::vector<QuadElement2D<FloatingDataType>*> neighbors_;
    LobattoTensorProductGrid2D<FloatingDataType> referenceGrid_;
    LobattoTensorProductGrid2D<FloatingDataType> physicalGrid_;
    LobattoTensorProductGrid2D<FloatingDataType> quadratureGrid_;
    std::map<std::string, Vector<FloatingDataType>> vecs_;
    std::map<std::string, Matrix<FloatingDataType>> mats_;
    FloatingDataType area_;

public:

    QuadElement2D(std::vector<SpaceVector2D<FloatingDataType>> points, int xOrder, int yOrder, int ID) :
        points_(points), boundaries_(4), neighbors_(4), referenceGrid_(xOrder, yOrder), physicalGrid_(xOrder, yOrder), quadratureGrid_(xOrder+1, yOrder+1), ID_(ID) {

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
            normals_[f] = normals_[f].normalize();
        }

        // Compute area
        SpaceVector2D<FloatingDataType> l0 = points_[1] - points_[0];
        SpaceVector2D<FloatingDataType> l1 = points_[3] - points_[0];
        area_ = -l0.y()*l1.x() + l0.x()*l1.y();

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
        // Matrix<FloatingDataType> dh_dxi = referenceGrid_.xPoly().derivative(quadratureGrid_.xPoints());
        // Matrix<FloatingDataType> h_eta = referenceGrid_.yPoly()(quadratureGrid_.yPoints());
        // Matrix<FloatingDataType> h_xi = referenceGrid_.xPoly()(quadratureGrid_.xPoints());
        // Matrix<FloatingDataType> dh_deta = referenceGrid_.yPoly().derivative(quadratureGrid_.yPoints());
        // Matrix<FloatingDataType> dpsi_dxi = kroneckerProduct(dh_dxi, h_eta);
        // Matrix<FloatingDataType> dpsi_deta = kroneckerProduct(h_xi, dh_deta);
        Matrix<FloatingDataType> dh_dxi = referenceGrid_.xPoly().derivative(referenceGrid_.xPoints());
        Matrix<FloatingDataType> h_eta = referenceGrid_.yPoly()(referenceGrid_.yPoints());
        Matrix<FloatingDataType> h_xi = referenceGrid_.xPoly()(referenceGrid_.xPoints());
        Matrix<FloatingDataType> dh_deta = referenceGrid_.yPoly().derivative(referenceGrid_.yPoints());
        Matrix<FloatingDataType> dpsi_dxi = kroneckerProduct(dh_dxi, h_eta);
        Matrix<FloatingDataType> dpsi_deta = kroneckerProduct(h_xi, dh_deta);

        // Compute metric terms (Algorithm 12.1)
        // Vector<FloatingDataType> dx_dxi(quadratureGrid_.size(), 0);
        // Vector<FloatingDataType> dx_deta(quadratureGrid_.size(), 0);
        // Vector<FloatingDataType> dy_dxi(quadratureGrid_.size(), 0);
        // Vector<FloatingDataType> dy_deta(quadratureGrid_.size(), 0);
        // Vector<FloatingDataType> dxi_dx(quadratureGrid_.size(), 0);
        // Vector<FloatingDataType> dxi_dy(quadratureGrid_.size(), 0);
        // Vector<FloatingDataType> deta_dx(quadratureGrid_.size(), 0);
        // Vector<FloatingDataType> deta_dy(quadratureGrid_.size(), 0);
        Vector<FloatingDataType> dx_dxi(referenceGrid_.size(), 0);
        Vector<FloatingDataType> dx_deta(referenceGrid_.size(), 0);
        Vector<FloatingDataType> dy_dxi(referenceGrid_.size(), 0);
        Vector<FloatingDataType> dy_deta(referenceGrid_.size(), 0);
        Vector<FloatingDataType> dxi_dx(referenceGrid_.size(), 0);
        Vector<FloatingDataType> dxi_dy(referenceGrid_.size(), 0);
        Vector<FloatingDataType> deta_dx(referenceGrid_.size(), 0);
        Vector<FloatingDataType> deta_dy(referenceGrid_.size(), 0);
        vecs_["jacobian"] = Vector<FloatingDataType>(referenceGrid_.size(), 0);
        for (auto k = 0; k < referenceGrid_.size(); k++) {
            for (auto i = 0; i < referenceGrid_.size(); i++) {
                SpaceVector2D<FloatingDataType> xi = referenceGrid_(i);
                SpaceVector2D<FloatingDataType> x = mapReference2Physical(xi);
                dx_dxi[k] += dpsi_dxi(i,k)*x.x();
                dx_deta[k] += dpsi_deta(i,k)*x.x();
                dy_dxi[k] += dpsi_dxi(i,k)*x.y();
                dy_deta[k] += dpsi_deta(i,k)*x.y();
            }
        }
        for (auto k = 0; k < referenceGrid_.size(); k++) {
            // vecs_["jacobian"][k] = vecs_["dx_dxi"][k]*vecs_["dy_deta"][k] - vecs_["dx_deta"][k]*vecs_["dy_dxi"][k];
            vecs_["jacobian"][k] = dx_deta[k]*dy_dxi[k] - dx_dxi[k]*dy_deta[k]; // ???
            // dxi_dx[k] = dy_deta[k] / vecs_["jacobian"][k];
            // dxi_dy[k] = -dx_deta[k] / vecs_["jacobian"][k];
            // deta_dx[k] = -dy_dxi[k] / vecs_["jacobian"][k];
            // deta_dy[k] = dx_dxi[k] / vecs_["jacobian"][k];
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
    std::vector<InternalBoundary<FloatingDataType>>& boundaries() { return boundaries_; }
    std::vector<QuadElement2D<FloatingDataType>*>& neighbors() { return neighbors_; }
    LobattoTensorProductGrid2D<FloatingDataType>& referenceGrid() { return referenceGrid_; }
    LobattoTensorProductGrid2D<FloatingDataType>& quadratureGrid() { return quadratureGrid_; }
    LobattoTensorProductGrid2D<FloatingDataType>& physicalGrid() { return physicalGrid_; }
    std::map<std::string, Vector<FloatingDataType>>& vectorMap() { return vecs_; }
    std::map<std::string, Matrix<FloatingDataType>>& matrixMap() { return mats_; }
    Vector<FloatingDataType>& vector(std::string name) { return vecs_[name]; }
    Matrix<FloatingDataType>& matrix(std::string name) { return mats_[name]; }
    int ID() { return ID_; }
    std::size_t size() { return referenceGrid_.size(); }
    FloatingDataType area() { return area_; }

    void setBoundary(int faceIndex, InternalBoundary<FloatingDataType>& boundary) {
        boundaries_[faceIndex] = &boundary;
    }

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

    Vector<FloatingDataType> computeVolumeIntegralInexact() {

        // Algorithm 16.9
        int N = referenceGrid_.size();
        Vector<FloatingDataType> R(N, 0);
        Vector<FloatingDataType>& flux_x = vecs_["fx"];
        Vector<FloatingDataType>& flux_y = vecs_["fy"];
        Matrix<FloatingDataType> W = referenceGrid_.basisWeights();

        // Compute grad(psi)
        Matrix<FloatingDataType> dpsi_x, dpsi_y;
        {
            Matrix<FloatingDataType> dh_dxi = referenceGrid_.xPoly().derivative(referenceGrid_.xPoints());
            Matrix<FloatingDataType> h_eta = referenceGrid_.yPoly()(referenceGrid_.yPoints());
            Matrix<FloatingDataType> h_xi = referenceGrid_.xPoly()(referenceGrid_.xPoints());
            Matrix<FloatingDataType> dh_deta = referenceGrid_.yPoly().derivative(referenceGrid_.yPoints());
            dpsi_y = kroneckerProduct(dh_dxi, h_eta);
            dpsi_x = kroneckerProduct(h_xi, dh_deta);
        }

        for (auto j = 0; j < N; j++) {
            for (auto i = 0; i < N; i++) {
                std::pair<int, int> ID = referenceGrid_.ID(j);
                FloatingDataType w_j = W(ID.first, ID.second);
                FloatingDataType J_j = vecs_["jacobian"][j];

                R[i] += w_j * J_j * (dpsi_x(i,j)*flux_x[j] + dpsi_y(i,j)*flux_y[j]);
            }
        }

        return R;

    }

    Vector<FloatingDataType> computeFluxIntegralInexact(int faceIndexInternal) {

        // Algorithm 16.10
        int N = referenceGrid_.size();
        Vector<FloatingDataType> R(N, 0);
        
        Vector<int> IS_internal = referenceGrid_.faceIndexSets()[faceIndexInternal];
        Vector<FloatingDataType>& q_internal = vecs_["q"];
        Vector<FloatingDataType>& flux_x_internal = vecs_["fx"];
        Vector<FloatingDataType>& flux_y_internal = vecs_["fy"];
        Vector<FloatingDataType>& u_x_internal = vecs_["ux"];
        Vector<FloatingDataType>& u_y_internal = vecs_["uy"];
        SpaceVector2D<FloatingDataType>& nhat = normals_[faceIndexInternal];
        Vector<FloatingDataType> w = (faceIndexInternal % 2) ? referenceGrid_.yWeights() : referenceGrid_.xWeights();
        Vector<FloatingDataType>& J = vecs_["jacobian"];

        QuadElement2D<FloatingDataType>& neighbor = *neighbors_[faceIndexInternal];
        int faceIndexExternal = neighborFaceIndex(faceIndexInternal);
        Vector<int> IS_external = neighbor.referenceGrid().faceIndexSets()[faceIndexExternal];
        Vector<FloatingDataType>& q_external = neighbor.vector("q");
        Vector<FloatingDataType>& flux_x_external = neighbor.vector("fx");
        Vector<FloatingDataType>& flux_y_external = neighbor.vector("fy");
        Vector<FloatingDataType>& u_x_external = neighbor.vector("ux");
        Vector<FloatingDataType>& u_y_external = neighbor.vector("uy");

        for (auto i = 0; i < IS_internal.size(); i++) {
            int ii = IS_internal[i];

            FloatingDataType u_internal = fabs(u_x_internal[i]*nhat.x() + u_y_internal[i]*nhat.y());
            FloatingDataType u_external = fabs(u_x_external[i]*nhat.x() + u_y_external[i]*nhat.y());
            FloatingDataType lambda = fmax(u_internal, u_external);
            FloatingDataType f_star_x = 0.5*(flux_x_internal[ii] + flux_x_external[ii] - lambda*(q_external[ii] - q_internal[ii])*nhat.x());
            FloatingDataType f_star_y = 0.5*(flux_y_internal[ii] + flux_y_external[ii] - lambda*(q_external[ii] - q_internal[ii])*nhat.y());

            R[ii] += w[i] * J[ii] * (f_star_x*nhat.x() + f_star_y*nhat.y());
        }

        return R;

    }

    Vector<FloatingDataType> computeVolumeIntegralExact() {

        // Algorithm 16.7
        int Q = quadratureGrid_.size();
        int N = referenceGrid_.size();
        Vector<FloatingDataType> R(N, 0);
        Vector<FloatingDataType>& flux_x = vecs_["fx"];
        Vector<FloatingDataType>& flux_y = vecs_["fy"];
        Matrix<FloatingDataType> psi = referenceGrid_.basisMatrix(quadratureGrid_);
        Matrix<FloatingDataType> W = quadratureGrid_.basisWeights();

        // Compute grad(psi)
        Matrix<FloatingDataType> dpsi_x, dpsi_y;
        {
            Matrix<FloatingDataType> dh_dxi = referenceGrid_.xPoly().derivative(quadratureGrid_.xPoints());
            Matrix<FloatingDataType> h_eta = referenceGrid_.yPoly()(quadratureGrid_.yPoints());
            Matrix<FloatingDataType> h_xi = referenceGrid_.xPoly()(quadratureGrid_.xPoints());
            Matrix<FloatingDataType> dh_deta = referenceGrid_.yPoly().derivative(quadratureGrid_.yPoints());
            dpsi_x = kroneckerProduct(dh_dxi, h_eta);
            dpsi_y = kroneckerProduct(h_xi, dh_deta);
        }

        for (auto k = 0; k < Q; k++) {
            Vector<FloatingDataType> f_x(Q, 0);
            Vector<FloatingDataType> f_y(Q, 0);

            for (auto j = 0; j < N; j++) {
                f_x[k] += psi(j,k) * flux_x[j];
                f_y[k] += psi(j,k) * flux_y[j];
            }

            for (auto i = 0; i < N; i++) {
                std::pair<int, int> ID = quadratureGrid_.ID(k);
                FloatingDataType w_k = W(ID.first, ID.second);
                FloatingDataType J_k = vecs_["jacobian"][k];
                R[i] += w_k * J_k * (dpsi_x(i,k)*f_x[k] + dpsi_y(i,k)*f_y[k]);
            }
        }

        return R;

    }

    Vector<FloatingDataType> computeFluxIntegralExact(int faceIndexInternal) {

        std::cout << "HERE" << std::endl;

        // Algorithm 16.8
        Vector<int> IS_referenceInternal = referenceGrid_.faceIndexSets()[faceIndexInternal];
        Vector<int> IS_quadratureInternal = quadratureGrid_.faceIndexSets()[faceIndexInternal];
        int N_l_internal = IS_referenceInternal.size();
        int Q_l_internal = IS_quadratureInternal.size();
        Vector<FloatingDataType>& q_internal = vecs_["q"];
        Vector<FloatingDataType> q_k_internal = q_internal.getFromIndexSet(IS_referenceInternal);
        Vector<FloatingDataType>& flux_x_internal = vecs_["fx"];
        Vector<FloatingDataType>& flux_y_internal = vecs_["fy"];

        QuadElement2D<FloatingDataType>& neighbor = *neighbors_[faceIndexInternal];
        int faceIndexExternal = neighborFaceIndex(faceIndexInternal);
        Vector<int> IS_referenceExternal = neighbor.referenceGrid().faceIndexSets()[faceIndexExternal];
        Vector<int> IS_quadratureExternal = neighbor.quadratureGrid().faceIndexSets()[faceIndexExternal];
        int N_l_external = IS_referenceExternal.size();
        int Q_l_external = IS_quadratureExternal.size();
        Vector<FloatingDataType>& q_external = neighbor.vector("q");
        Vector<FloatingDataType> q_k_external = q_external.getFromIndexSet(IS_referenceExternal);
        Vector<FloatingDataType>& flux_x_external = neighbor.vector("fx");
        Vector<FloatingDataType>& flux_y_external = neighbor.vector("fy");

        Vector<FloatingDataType> w_l = (faceIndexInternal % 2) ? quadratureGrid_.yWeights() : quadratureGrid_.xWeights();
        Vector<FloatingDataType> J_l = vecs_["jacobian"].getFromIndexSet(IS_quadratureInternal);
        Matrix<FloatingDataType> psi = referenceGrid_.basisMatrix(quadratureGrid_);
        SpaceVector2D<FloatingDataType>& nhat = normals_[faceIndexInternal];

        Vector<FloatingDataType> R(N_l_internal, 0);

        for (auto k = 0; k < Q_l_internal; k++) {
            int kk = IS_quadratureInternal[k];

            Vector<FloatingDataType> fL_x(Q_l_internal, 0);
            Vector<FloatingDataType> fL_y(Q_l_internal, 0);
            Vector<FloatingDataType> fR_x(Q_l_internal, 0);
            Vector<FloatingDataType> fR_y(Q_l_internal, 0);

            for (auto j = 0; j < N_l_internal; j++) {
                int jj = IS_referenceInternal[j];

                fL_x[kk] += psi(jj,kk) * flux_x_internal[jj];
                fL_y[kk] += psi(jj,kk) * flux_y_internal[jj];
                fR_x[kk] += psi(jj,kk) * flux_x_external[jj];
                fR_y[kk] += psi(jj,kk) * flux_y_external[jj];
                
            }
            
            FloatingDataType lambda = 1.0; // TODO: What is this??
            FloatingDataType f_star_x = 0.5*(fL_x[kk] + fR_x[kk] - lambda*(q_k_external - q_k_internal)*nhat.x());

            for (auto i = 0; i < N_l_internal; i++) {
                int ii = IS_referenceInternal[i];

                // R[i] -= w_l[k] * J_l[k] * psi(ii,kk) * (nhat.x()*f_star_x)
            }

        }

        return R;

    }

};

} // NAMESPACE : HydroForest

#endif // ELEMENT_2D_HPP_