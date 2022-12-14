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
class CGDerivativeMatrix : public Matrix<NumericalType> {

public:

    CGDerivativeMatrix(LagrangePolynomial& basisFunctionPoly, Grid1DBase<NumericalType>& nodalPointGrid, Grid1DBase<NumericalType>& quadraturePointGrid) :
        Matrix<NumericalType>(basisFunctionPoly.order, basisFunctionPoly.order, 0) {

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

template<typename NumericalType>
class CGLaplacianMatrix : public Matrix<NumericalType> {

public:

    CGLaplacianMatrix(LagrangePolynomial& basisFunctionPoly, Grid1DBase<NumericalType>& nodalPointGrid, Grid1DBase<NumericalType>& quadraturePointGrid) :
        Matrix<NumericalType>(basisFunctionPoly.order, basisFunctionPoly.order, 0) {

        Matrix<NumericalType> dL_ik = basisFunctionPoly.derivative(quadraturePointGrid.getPoints());
        Matrix<NumericalType> dL_jk = basisFunctionPoly.derivative(quadraturePointGrid.getPoints());

        for (auto k = 0; k < quadraturePointGrid.getNPoints(); k++) {
            for (auto j = 0; j < nodalPointGrid.getNPoints(); j++) {
                for (auto i = 0; i < nodalPointGrid.getNPoints(); i++) {
                    NumericalType w_k = quadraturePointGrid.getWeights()[k];
                    NumericalType dphi_ik = dL_ik(i,k);
                    NumericalType dphi_jk = dL_jk(j,k);
                    this->operator()(i,j) += w_k * dphi_ik * dphi_jk;
                }
            }
        }

    }

};

class CGIDMatrix : public Matrix<int> {

public:

    CGIDMatrix(std::size_t nOrder, std::size_t nElements) :
        Matrix<int>(nOrder + 1, nElements, 0) {

        int ID_ij = 0;
        for (auto e = 0; e < nElements; e++) {
            for (auto i = 0; i < nOrder+1; i++) {
                this->operator()(i, e) = ID_ij++;
            }
            ID_ij--;
        }

    }

};

template<typename NumericalType>
struct CGDirectStiffnessSummationOperator {

    std::vector<Element1D<NumericalType>>& elements;
    Matrix<int>& IDMatrix;

    CGDirectStiffnessSummationOperator(std::vector<Element1D<NumericalType>>& elements, Matrix<int>& IDMatrix) :
        elements(elements), IDMatrix(IDMatrix)
            {}

    Matrix<NumericalType> operate(Matrix<NumericalType>& elementMatrix) {
        int nElements = (int) elements.size();
        int nOrder = (int) IDMatrix.nRows() - 1;
        int nPoints = nElements * nOrder + 1;

        Matrix<NumericalType> M(nPoints, nPoints, 0);

        for (auto e = 0; e < nElements; e++) {
            for (auto j = 0; j < nOrder + 1; j++) {
                int J = IDMatrix(j, e);
                for (auto i = 0; i < nOrder + 1; i++) {
                    int I = IDMatrix(i, e);
                    M(I,J) += elementMatrix(i,j);
                }
            }
        }

        return M;
    }

    Vector<NumericalType> operate(Vector<NumericalType>& elementVector) {
        int nElements = (int) elements.size();
        int nOrder = (int) IDMatrix.nRows() - 1;
        int nPoints = nElements * nOrder + 1;

        Vector<NumericalType> r(nPoints, 0);

        for (auto e = 0; e < nElements; e++) {
            for (auto i = 0; i < nOrder; i++) {
                int I = IDMatrix(i,e);
                r[I] += elementVector[i];
            }
        }

        return r;
    }

    Matrix<NumericalType> operateWithMetricTerm(Matrix<NumericalType>& elementMatrix) {
        int nElements = (int) elements.size();
        int nOrder = (int) IDMatrix.nRows() - 1;
        int nPoints = nElements * nOrder + 1;

        Matrix<NumericalType> M(nPoints, nPoints, 0);

        for (auto e = 0; e < nElements; e++) {
            double dx = elements[e].xUpper() - elements[e].xLower();
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

    Matrix<NumericalType> operateWithInverseMetricTerm(Matrix<NumericalType>& elementMatrix) {
        int nElements = (int) elements.size();
        int nOrder = (int) IDMatrix.nRows() - 1;
        int nPoints = nElements * nOrder + 1;

        Matrix<NumericalType> M(nPoints, nPoints, 0);

        for (auto e = 0; e < nElements; e++) {
            double dx = elements[e].xUpper() - elements[e].xLower();
            for (auto j = 0; j < nOrder + 1; j++) {
                int J = IDMatrix(j, e);
                for (auto i = 0; i < nOrder + 1; i++) {
                    int I = IDMatrix(i, e);
                    M(I,J) += (2.0/dx)*elementMatrix(i,j);
                }
            }
        }

        return M;
    }

    Vector<NumericalType> operateWithMetricTerm(Vector<NumericalType>& elementVector) {
        int nElements = (int) elements.size();
        int nOrder = (int) IDMatrix.nRows() - 1;
        int nPoints = nElements * nOrder + 1;

        Vector<NumericalType> r(nPoints, 0);

        for (auto e = 0; e < nElements; e++) {
            double dx = elements[e].xUpper() - elements[e].xLower();
            for (auto i = 0; i < nOrder; i++) {
                int I = IDMatrix(i,e);
                r[I] += (dx/2.0)*elementVector[i];
            }
        }

        return r;
    }

    Matrix<NumericalType> operateWithString(std::string matrixName) {
        int nElements = (int) elements.size();
        int nOrder = (int) IDMatrix.nRows() - 1;
        int nPoints = nElements * nOrder + 1;

        Matrix<NumericalType> M(nPoints, nPoints, 0);

        for (auto e = 0; e < nElements; e++) {
            for (auto j = 0; j < nOrder + 1; j++) {
                int J = IDMatrix(j, e);
                for (auto i = 0; i < nOrder + 1; i++) {
                    int I = IDMatrix(i, e);
                    M(I,J) += elements[e].matrices()[matrixName](i,j);
                }
            }
        }

        return M;
    }

    Vector<NumericalType> operateWithStringAndMetricTermVector(std::string vectorName) {
        int nElements = (int) elements.size();
        int nOrder = (int) IDMatrix.nRows() - 1;
        int nPoints = nElements * nOrder + 1;

        Vector<NumericalType> r(nPoints, 0);

        for (auto e = 0; e < nElements; e++) {
            for (auto i = 0; i < nOrder; i++) {
                int I = IDMatrix(i,e);
                r[I] += elements[e].vectors()[vectorName][i];
            }
        }

        return r;
    }

    Matrix<NumericalType> operateWithStringAndMetricTermMatrix(std::string matrixName) {
        int nElements = (int) elements.size();
        int nOrder = (int) IDMatrix.nRows() - 1;
        int nPoints = nElements * nOrder + 1;

        Matrix<NumericalType> M(nPoints, nPoints, 0);

        for (auto e = 0; e < nElements; e++) {
            double dx = elements[e].xUpper() - elements[e].xLower();
            for (auto j = 0; j < nOrder + 1; j++) {
                int J = IDMatrix(j, e);
                for (auto i = 0; i < nOrder + 1; i++) {
                    int I = IDMatrix(i, e);
                    M(I,J) += (dx/2.0)*elements[e].matrices()[matrixName](i,j);
                }
            }
        }

        return M;
    }

    Vector<NumericalType> operateWithStringAndMetricTerm(std::string vectorName) {
        int nElements = (int) elements.size();
        int nOrder = (int) IDMatrix.nRows() - 1;
        int nPoints = nElements * nOrder + 1;

        Vector<NumericalType> r(nPoints, 0);

        for (auto e = 0; e < nElements; e++) {
            double dx = elements[e].xUpper() - elements[e].xLower();
            for (auto i = 0; i < nOrder; i++) {
                int I = IDMatrix(i,e);
                r[I] += (dx/2.0)*elements[e].vectors()[vectorName][i];
            }
        }

        return r;
    }

};

} // NAMESPACE : HydroForest

#endif // CG_FINITE_ELEMENTS_HPP_