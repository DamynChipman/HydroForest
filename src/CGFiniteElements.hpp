#ifndef CG_FINITE_ELEMENTS_HPP_
#define CG_FINITE_ELEMENTS_HPP_

#include "HydroForestApp.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "Polynomial.hpp"
#include "Element1D.hpp"

namespace HydroForest {

template<typename NumericalType>
class CGMassMatrix : public Matrix<NumericalType> {

public:

    CGMassMatrix(LagrangePolynomial& basisFunctionPoly, Grid1DBase<NumericalType>& nodalPointGrid, Grid1DBase<NumericalType>& quadraturePointGrid) :
        Matrix<NumericalType>(basisFunctionPoly.order, basisFunctionPoly.order) {

        for (auto& i : this->data_) {
            i = 0;
        }

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
class CGDerivativeMatrix : public Matrix<NumericalType> {

public:

    CGDerivativeMatrix(LagrangePolynomial& basisFunctionPoly, Grid1DBase<NumericalType>& nodalPointGrid, Grid1DBase<NumericalType>& quadraturePointGrid) :
        Matrix<NumericalType>(basisFunctionPoly.order, basisFunctionPoly.order) {

        for (auto& i : this->data_) i = 0;

        Matrix<NumericalType> L_ik = basisFunctionPoly(quadraturePointGrid.getPoints());
        Matrix<NumericalType> dL_jk = basisFunctionPoly.derivative(quadraturePointGrid.getPoints());

        for (auto k = 0; k < quadraturePointGrid.getNPoints(); k++) {
            for (auto j = 0; j < nodalPointGrid.getNPoints(); j++) {
                for (auto i = 0; i < nodalPointGrid.getNPoints(); i++) {
                    NumericalType w_k = quadraturePointGrid.getWeights()[k];
                    NumericalType phi_ik = L_ik(i,k);
                    NumericalType dphi_jk = dL_jk(j,k);
                    this->operator()(i, j) += w_k * phi_ik * dphi_jk;
                }
            }
        }

    }

};

class CGIDMatrix : public Matrix<int> {

public:

    CGIDMatrix(std::size_t nOrder, std::size_t nElements) :
        Matrix<int>(nOrder + 1, nElements) {

        for (auto& i : this->data_) i = 0;

        int ID_ij = 0;
        for (auto j = 0; j < nElements; j++) {
            for (auto i = 0; i < nOrder+1; i++) {
                this->operator()(i, j) = ID_ij++;
            }
            ID_ij--;
        }

    }

};

template<typename NumericalType>
Matrix<NumericalType> directStiffnessSummation(std::vector<Element1D<NumericalType>>& elements, Matrix<int>& IDMatrix, Matrix<NumericalType>& elementMatrix) {
    
    int nElements = (int) elements.size();
    int nOrder = (int) IDMatrix.nRows() - 1;
    int nPoints = nElements * nOrder + 1;

    Matrix<NumericalType> M(nPoints, nPoints, 0);

    for (auto e = 0; e < nElements; e++) {
        double dx = elements[e].grid().getUpperBound() - elements[e].grid().getLowerBound();
        for (auto j = 0; j < nOrder + 1; j++) {
            int J = IDMatrix(j, e);
            for (auto i = 0; i < nOrder + 1; i++) {
                int I = IDMatrix(i, e);
                M(I,J) += (dx/2.0)*elementMatrix(i,j);
            }
        }
    }

    return M;

}

} // NAMESPACE : HydroForest

#endif // CG_FINITE_ELEMENTS_HPP_