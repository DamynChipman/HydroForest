#ifndef DG_FINITE_ELEMENTS_HPP_
#define DG_FINITE_ELEMENTS_HPP_

#include "HydroForestApp.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "Polynomial.hpp"
#include "Element1D.hpp"
#include "Element2D.hpp"
#include "Mesh1D.hpp"
#include "Grid2D.hpp"

namespace HydroForest {

template<typename NumericalType>
class DGMassMatrix : public Matrix<NumericalType> {

public:

    DGMassMatrix(LagrangePolynomial& basisFunctionPoly, Grid1DBase<NumericalType>& nodalPointGrid, Grid1DBase<NumericalType>& quadraturePointGrid) :
        Matrix<NumericalType>(basisFunctionPoly.order, basisFunctionPoly.order, 0) {

        Matrix<NumericalType> L_ik = basisFunctionPoly(quadraturePointGrid.getPoints());
        Matrix<NumericalType> L_jk = basisFunctionPoly(quadraturePointGrid.getPoints());

        for (auto k = 0; k < quadraturePointGrid.getNPoints(); k++) {
            for (auto j = 0; j < nodalPointGrid.getNPoints(); j++) {
                for (auto i = 0; i < nodalPointGrid.getNPoints(); i++) {
                    NumericalType w_k = quadraturePointGrid.getWeights()[k];
                    NumericalType phi_ik = L_ik(i,k);
                    NumericalType phi_jk = L_jk(j,k);
                    this->operator()(i, j) += w_k * phi_ik * phi_jk;
                }
            }
        }

    }

};

template<typename NumericalType>
class DGMassMatrix2D : public Matrix<NumericalType> {

public:

    DGMassMatrix2D(LobattoTensorProductGrid2D<NumericalType>& nodalPointGrid, LobattoTensorProductGrid2D<NumericalType>& quadraturePointGrid) :
        Matrix<NumericalType>(nodalPointGrid.xSize()*nodalPointGrid.ySize(), nodalPointGrid.xSize()*nodalPointGrid.ySize(), 0) {

        //
        int M = quadraturePointGrid.xSize()*quadraturePointGrid.ySize();
        int N = nodalPointGrid.xSize()*nodalPointGrid.ySize();
        Matrix<NumericalType> phi = nodalPointGrid.basisMatrix(quadraturePointGrid);
        Matrix<NumericalType> W = quadraturePointGrid.basisWeights();

        for (auto k = 0; k < M; k++) {
            for (auto i = 0; i < N; i++) {
                for (auto j = 0; j < N; j++) {
                    std::pair<int, int> ID = quadraturePointGrid.ID(k);
                    NumericalType w_k = W(ID.first, ID.second);
                    this->operator()(i, j) += w_k * phi(i,k) * phi(j,k);
                }
            }
        }

    }

    DGMassMatrix2D(QuadElement2D<NumericalType>& element) :
        Matrix<NumericalType>(element.referenceGrid().xSize()*element.referenceGrid().ySize(), element.referenceGrid().xSize()*element.referenceGrid().ySize(), 0) {

        // Algorithm 12.6
        LobattoTensorProductGrid2D<NumericalType>& nodalPointGrid = element.referenceGrid();
        LobattoTensorProductGrid2D<NumericalType>& quadraturePointGrid = element.quadratureGrid();
        int M = quadraturePointGrid.xSize()*quadraturePointGrid.ySize();
        int N = nodalPointGrid.xSize()*nodalPointGrid.ySize();
        Matrix<NumericalType> phi = nodalPointGrid.basisMatrix(quadraturePointGrid);
        Matrix<NumericalType> W = quadraturePointGrid.basisWeights();

        for (auto k = 0; k < M; k++) {
            for (auto i = 0; i < N; i++) {
                for (auto j = 0; j < N; j++) {
                    std::pair<int, int> ID = quadraturePointGrid.ID(k);
                    NumericalType w_k = W(ID.first, ID.second);
                    NumericalType J_k = element.vector("jacobian")[k];
                    this->operator()(i, j) += w_k * J_k * phi(i,k) * phi(j,k);
                }
            }
        }

    }

};

template<typename NumericalType>
class DGDerivativeMatrix : public Matrix<NumericalType> {

public:

    DGDerivativeMatrix(LagrangePolynomial& basisFunctionPoly, Grid1DBase<NumericalType>& nodalPointGrid, Grid1DBase<NumericalType>& quadraturePointGrid) :
        Matrix<NumericalType>(basisFunctionPoly.order, basisFunctionPoly.order, 0) {

        // for (auto& i : this->data_) i = 0;

        Matrix<NumericalType> dL_ik = basisFunctionPoly.derivative(quadraturePointGrid.getPoints());
        Matrix<NumericalType> L_jk = basisFunctionPoly(quadraturePointGrid.getPoints());

        for (auto k = 0; k < quadraturePointGrid.getNPoints(); k++) {
            for (auto j = 0; j < nodalPointGrid.getNPoints(); j++) {
                for (auto i = 0; i < nodalPointGrid.getNPoints(); i++) {
                    NumericalType w_k = quadraturePointGrid.getWeights()[k];
                    NumericalType dphi_ik = dL_ik(i,k);
                    NumericalType phi_jk = L_jk(j,k);
                    this->operator()(i, j) += w_k * dphi_ik * phi_jk;
                }
            }
        }

    }

};

template<typename NumericalType>
class DGDerivativeMatrix2D : public Matrix<NumericalType> {

public:

    // DGDerivativeMatrix2D(LobattoTensorProductGrid2D<NumericalType>& nodalPointGrid, LobattoTensorProductGrid2D<NumericalType>& quadraturePointGrid) :
    //     Matrix<NumericalType>(nodalPointGrid.xSize()*nodalPointGrid.ySize(), nodalPointGrid.xSize()*nodalPointGrid.ySize(), 0) {

    //     //
    //     int M = quadraturePointGrid.xSize()*quadraturePointGrid.ySize();
    //     int N = nodalPointGrid.xSize()*nodalPointGrid.ySize();
    //     Matrix<NumericalType> phi = nodalPointGrid.basisMatrix(quadraturePointGrid);
    //     Matrix<NumericalType> dphi = nodalPointGrid.basisDerivativeMatrix(quadraturePointGrid);
    //     Matrix<NumericalType> W = quadraturePointGrid.basisWeights();

    //     for (auto k = 0; k < M; k++) {
    //         for (auto i = 0; i < N; i++) {
    //             for (auto j = 0; j < N; j++) {
    //                 std::pair<int, int> ID = quadraturePointGrid.ID(k);
    //                 NumericalType w_k = W(ID.first, ID.second);
    //                 this->operator()(i, j) += w_k * dphi(i,k) * phi(j,k);
    //             }
    //         }
    //     }

    // }

    DGDerivativeMatrix2D(QuadElement2D<NumericalType>& element, int dim) :
        Matrix<NumericalType>(element.referenceGrid().xSize()*element.referenceGrid().ySize(), element.referenceGrid().xSize()*element.referenceGrid().ySize(), 0) {

        // Algorithm 16.1
        LobattoTensorProductGrid2D<NumericalType>& nodalPointGrid = element.referenceGrid();
        LobattoTensorProductGrid2D<NumericalType>& quadraturePointGrid = element.quadratureGrid();
        int Q = quadraturePointGrid.xSize()*quadraturePointGrid.ySize();
        int N = nodalPointGrid.xSize()*nodalPointGrid.ySize();
        Matrix<NumericalType> psi = nodalPointGrid.basisMatrix(quadraturePointGrid);
        Matrix<NumericalType> W = quadraturePointGrid.basisWeights();

        Matrix<NumericalType> dpsi;
        if (dim == 0) {
            Matrix<NumericalType> dh_dxi = nodalPointGrid.xPoly().derivative(quadraturePointGrid.xPoints());
            Matrix<NumericalType> h_eta = nodalPointGrid.yPoly()(quadraturePointGrid.yPoints());
            dpsi = kroneckerProduct(dh_dxi, h_eta);
        }
        else if (dim == 1) {
            Matrix<NumericalType> h_xi = nodalPointGrid.xPoly()(quadraturePointGrid.xPoints());
            Matrix<NumericalType> dh_deta = nodalPointGrid.yPoly().derivative(quadraturePointGrid.yPoints());
            dpsi = kroneckerProduct(h_xi, dh_deta);
        }

        for (auto k = 0; k < Q; k++) {
            for (auto i = 0; i < N; i++) {
                for (auto j = 0; j < N; j++) {
                    std::pair<int, int> ID = quadraturePointGrid.ID(k);
                    NumericalType w_k = W(ID.first, ID.second);
                    NumericalType J_k = element.vector("jacobian")[k];
                    this->operator()(i, j) += w_k * J_k * dpsi(i,k) * psi(j,k);
                }
            }
        }

    }

};

class DGIDMatrix : public Matrix<int> {

public:

    DGIDMatrix(std::size_t nOrder, std::size_t nElements) :
        Matrix<int>(nOrder + 1, nElements, 0) {

        int ID_ij = 0;
        for (auto e = 0; e < nElements; e++) {
            for (auto i = 0; i < nOrder+1; i++) {
                this->operator()(i,e) = ID_ij++;
            }
        }

    }

};

template<typename NumericalType>
class DGFluxMatrix : public Matrix<NumericalType> {

public:

    DGFluxMatrix(std::size_t order) :
        Matrix<NumericalType>(order+1, order+1, 0) {

        this->operator()(0, 0) = -1;
        this->operator()(order, order) = 1;

    }

};

template<typename NumericalType>
class DGFluxMatrix2D : public Matrix<NumericalType> {

public:

    DGFluxMatrix2D(QuadElement2D<NumericalType>& element, int dim, int faceIndex) :
        Matrix<NumericalType>(element.referenceGrid().xSize()*element.referenceGrid().ySize(), element.referenceGrid().xSize()*element.referenceGrid().ySize(), 0) {

        // Algorithm 16.4
        Vector<int> referenceGridFaceIndexSet = element.referenceGrid().faceIndexSets()[faceIndex];
        Vector<int> quadratureGridFaceIndexSet = element.quadratureGrid().faceIndexSets()[faceIndex];
        int N_l = referenceGridFaceIndexSet.size();
        int Q_l = quadratureGridFaceIndexSet.size();
        Vector<NumericalType> w_l = (faceIndex % 2) ? element.quadratureGrid().yWeights() : element.quadratureGrid().xWeights();
        Vector<NumericalType> J_l = element.vector("jacobian").getFromIndexSet(quadratureGridFaceIndexSet);
        Matrix<NumericalType> psi = element.referenceGrid().basisMatrix(element.quadratureGrid());
        NumericalType nhat = element.normals()[faceIndex][dim];

        // int N_l;
        // int Q_l;
        // Vector<NumericalType> w_l;
        // Vector<int> jacobianIndexSet;
        // if (faceIndex == 0) { // South
        //     N_l = element.referenceGrid().xSize();
        //     Q_l = element.quadratureGrid().xSize();
        //     w_l = element.quadratureGrid().xWeights();
        //     jacobianIndexSet = Vector<int>(Q_l);
        //     for (auto i = 0; i < Q_l; i++) jacobianIndexSet[i] = i;
        // }
        // else if (faceIndex == 1) { // East
        //     N_l = element.referenceGrid().ySize();
        //     Q_l = element.quadratureGrid().ySize();
        //     w_l = element.quadratureGrid().yWeights();
        //     jacobianIndexSet = Vector<int>(Q_l);
        //     for (auto j = 0; j < Q_l; j++) jacobianIndexSet[j] = (j+1)*(element.quadratureGrid().xSize()) - 1;
        // }
        // else if (faceIndex == 2) { // North
        //     N_l = element.referenceGrid().xSize();
        //     Q_l = element.quadratureGrid().xSize();
        //     w_l = element.quadratureGrid().xWeights();
        //     jacobianIndexSet = Vector<int>(Q_l);
        //     for (auto i = 0; i < Q_l; i++) jacobianIndexSet[i] = (element.quadratureGrid().ySize()-1)*element.quadratureGrid().ySize() + i;
        // }
        // else if (faceIndex == 3) { // West
        //     N_l = element.referenceGrid().ySize();
        //     Q_l = element.quadratureGrid().ySize();
        //     w_l = element.quadratureGrid().yWeights();
        //     jacobianIndexSet = Vector<int>(Q_l);
        //     for (auto j = 0; j < Q_l; j++) jacobianIndexSet[j] = (j)*(element.quadratureGrid().xSize());
        // }
        // else {
        //     std::string errorMessage = "[HydroForest::DGFluxMatrix2D] Invalid `faceIndex`. Options are [0,3]:\n";
        //     errorMessage += "\tfaceIndex = " + std::to_string(faceIndex);
        //     throw std::invalid_argument(errorMessage);
        // }
        // Vector<NumericalType> J_l = element.vector("jacobian").getFromIndexSet(jacobianIndexSet);        

        for (auto k = 0; k < Q_l; k++) {
            for (auto i = 0; i < N_l; i++) {
                for (auto j = 0; j < N_l; j++) {
                    int kk = quadratureGridFaceIndexSet[k];
                    int ii = referenceGridFaceIndexSet[i];
                    int jj = referenceGridFaceIndexSet[j];
                    this->operator()(ii,jj) += w_l[k] * J_l[k] * psi(i,k) * psi(j,k) * nhat; // WORKING HERE
                }
            }
        }

    }

};

template<typename NumericalType>
class DGGlobalCenteredFluxMatrix : public Matrix<NumericalType> {

public:

    DGGlobalCenteredFluxMatrix(std::vector<Element1D<NumericalType>>& elements, Matrix<int>& IDMatrix, BoundaryConditionType leftBoundaryType, BoundaryConditionType rightBoundaryType) :
        Matrix<NumericalType>(elements.size()*(IDMatrix.nRows()), elements.size()*(IDMatrix.nRows()), 0) {

        //
        int N = IDMatrix.nRows()-1;
        for (auto e = 0; e < elements.size(); e++) {
            int L, R, I, J;

            // Identify left and right elements
            L = e-1;
            R = e+1;
            if (e == 0) {
                if (leftBoundaryType == BoundaryConditionType::Periodic) L = elements.size()-1;
                else L = -1;
            }
            else if (e == elements.size()-1) {
                if (rightBoundaryType == BoundaryConditionType::Periodic) R = 0;
                else R = -1;
            }

            if (L != -1) {
                I = IDMatrix(0, e);
                J = IDMatrix(N, L);
                this->operator()(I,I) = -0.5;
                this->operator()(I,J) = -0.5;
            }
            else {
                I = 0;
                this->operator()(I,I) = 1.0;
            }

            if (R != -1) {
                I = IDMatrix(N, e);
                J = IDMatrix(0, R);
                this->operator()(I,I) = 0.5;
                this->operator()(I,J) = 0.5;
            }
            else {
                I = this->nRows()-1;
                this->operator()(I,I) = 1.0;
            }
        }

    }

};

template<typename NumericalType>
class DGGlobalRusanovFluxVector : public Vector<NumericalType> {

public:

    DGGlobalRusanovFluxVector(std::vector<Element1D<NumericalType>>& elements, Matrix<int>& IDMatrix, Vector<NumericalType>& q_global, Vector<NumericalType>& f_global, NumericalType lambda) :
        Vector<NumericalType>(elements.size()*(IDMatrix.nRows()), 0) {

        //
        int N = IDMatrix.nRows()-1;
        double schemeFlag = 1.0;
        for (auto e = 0; e < elements.size(); e++) {
            int L = e;
            int R = (e + 1) % elements.size();
            int I = IDMatrix(N, L);
            int J = IDMatrix(0, R);
            NumericalType f_star = 0.5*(f_global[I] + f_global[J] - schemeFlag*lambda*(q_global[J] - q_global[I]));
            this->data_[I] = f_star;
            this->data_[J] = f_star;
        }

    }

};

template<typename NumericalType>
class DGGlobalFluxVector : public Vector<NumericalType> {

protected:

    double schemeFlag_;

public:

    // For Rusanov, supply f_global
    DGGlobalFluxVector(std::vector<Element1D<NumericalType>>& elements, Matrix<int>& IDMatrix, Vector<NumericalType>& q_global, Vector<NumericalType>& f_global, NumericalType lambda, double scheme=1.0) :
        Vector<NumericalType>(elements.size()*(IDMatrix.nRows()), 0), schemeFlag_(scheme) {

        //
        int N = IDMatrix.nRows()-1;
        for (auto e = 0; e < elements.size(); e++) {
            int L = e;
            int R = (e + 1) % elements.size();
            int I = IDMatrix(N, L);
            int J = IDMatrix(0, R);
            NumericalType f_star = 0.5*(f_global[I] + f_global[J] - schemeFlag_*lambda*(q_global[J] - q_global[I]));
            this->data_[I] = f_star;
            this->data_[J] = f_star;
        }

    }

    // For Centered, only supply q_global
    DGGlobalFluxVector(std::vector<Element1D<NumericalType>>& elements, Matrix<int>& IDMatrix, Vector<NumericalType>& q_global) :
        Vector<NumericalType>(elements.size()*(IDMatrix.nRows()), 0), schemeFlag_(0.0) {

        //
        int N = IDMatrix.nRows()-1;
        for (auto e = 0; e < elements.size(); e++) {
            int L = e;
            int R = (e + 1) % elements.size();
            int I = IDMatrix(N, L);
            int J = IDMatrix(0, R);
            NumericalType f_star = 0.5*(q_global[J] - q_global[I]);
            this->data_[I] = f_star;
            this->data_[J] = f_star;
        }

    }

};

template<typename NumericalType>
struct DGDirectStiffnessSummationOperator {

    std::vector<Element1D<NumericalType>>& elements;
    Matrix<int>& IDMatrix;

    DGDirectStiffnessSummationOperator(std::vector<Element1D<NumericalType>>& elements, Matrix<int>& IDMatrix) :
        elements(elements), IDMatrix(IDMatrix)
            {}

    Matrix<NumericalType> operate(Matrix<NumericalType>& elementMatrix) {
        int nElements = (int) elements.size();
        std::vector<Matrix<NumericalType>> diag(nElements);
        for (auto e = 0; e < nElements; e++) {
            diag[e] = elementMatrix;
        }
        return blockDiagonalMatrix(diag);
    }

    Matrix<NumericalType> operateWithMetricTerm(Matrix<NumericalType>& elementMatrix) {
        int nElements = (int) elements.size();
        std::vector<Matrix<NumericalType>> diag(nElements);
        for (auto e = 0; e < nElements; e++) {
            diag[e] = elementMatrix;
            double dx = elements[e].xUpper() - elements[e].xLower();
            diag[e] *= (dx/2.0);
        }
        return blockDiagonalMatrix(diag);
    }

};

} // NAMESPACE : HydroForest

#endif // DG_FINITE_ELEMENTS_HPP_