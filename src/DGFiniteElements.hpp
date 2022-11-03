#ifndef DG_FINITE_ELEMENTS_HPP_
#define DG_FINITE_ELEMENTS_HPP_

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