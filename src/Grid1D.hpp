#ifndef GRID_1D_HPP_
#define GRID_1D_HPP_

#include <vector>

#include "Polynomial.hpp"
#include "VTKBuilder.hpp"
#include "NewtonRaphsonSolver.hpp"

namespace HydroForest {

template<typename FloatingDataType>
class Grid1DBase {

public:

    virtual std::vector<FloatingDataType> getPoints() = 0;
    virtual std::size_t getNPoints() = 0;
    virtual FloatingDataType getLowerBound() = 0;
    virtual FloatingDataType getUpperBound() = 0;

};

template<typename FloatingDataType>
class UniformGrid1D : public Grid1DBase<FloatingDataType>, public RectilinearGridNodeBase, public DataArrayNodeBase {

public:

    UniformGrid1D(FloatingDataType lowerBound, FloatingDataType upperBound, std::size_t nPoints) : lowerBound_(lowerBound), upperBound_(upperBound), nPoints_(nPoints), points_(nPoints) {
        points_.resize(nPoints_);
        gridSpacing_ = (upperBound_ - lowerBound_) / (nPoints_ - 1);
        for (auto i = 0; i < nPoints_; i++) {
            points_[i] = lowerBound_ + i * gridSpacing_;
        }
    }

    std::vector<FloatingDataType> getPoints() { return points_; }
    std::size_t getNPoints() { return points_.size(); }
    FloatingDataType getLowerBound() { return lowerBound_; }
    FloatingDataType getUpperBound() { return upperBound_; }

    std::string getWholeExtent() {
        return "0 " + std::to_string(nPoints_) + " 0 0 0 0";
    }

    std::string getExtent() {
        return "0 " + std::to_string(nPoints_) + " 0 0 0 0";
    }

    std::string getType() {
        return "Float32";
    }

    std::string getName() {
        return "UniformGrid1d";
    }

    std::string getNumberOfComponents() {
        return "1";
    }

    std::string getFormat() {
        return "ascii";
    }

    std::string getRangeMin() {
        return std::to_string(lowerBound_);
    }

    std::string getRangeMax() {
        return std::to_string(upperBound_);
    }

    std::string getData() {
        if (nPoints_ == 0) {
            return "0";
        }
        else {
            std::string pointsAsString = "";
            for (auto& p : points_) pointsAsString += std::to_string(p) + " ";
            return pointsAsString;
        }
    }

private:

    FloatingDataType lowerBound_;
    FloatingDataType upperBound_;
    FloatingDataType gridSpacing_;
    std::size_t nPoints_;
    std::vector<FloatingDataType> points_;

};

template <typename FloatingDataType>
class ChebyshevGrid1D : public Grid1DBase<FloatingDataType> {

public:

    ChebyshevGrid1D(std::size_t order) : order_(order), points_(order+1), weights_(order+1) {

        for (int i = 0; i <= order_; i++) {
            points_[i] = cos((2.0*i + 1.0)/(2.0*order_ + 2.0) * M_PI);
        }

        for (int i = 0; i<= order_; i++) {
            weights_[i] = (M_PI) / (order_ + 1.0);
        }

    }

    double weightFunction(double x) {
        return (1.0) / (sqrt(1.0 - pow(x, 2)));
    }

    std::vector<FloatingDataType> getPoints() { return points_; }
    std::size_t getNPoints() { return points_.size(); }
    FloatingDataType getLowerBound() { return -1.0; }
    FloatingDataType getUpperBound() { return 1.0; }

private:

    std::size_t order_;
    std::vector<FloatingDataType> points_;
    std::vector<FloatingDataType> weights_;

};

template <typename FloatingDataType>
class LegendreGrid1D : public Grid1DBase<FloatingDataType> {

public:

    LegendrePolynomial poly;

    struct RootSolver : public NewtonRaphsonSolver<LegendreGrid1D> {

        LegendreGrid1D& actualClass;

        RootSolver(LegendreGrid1D& actualClass, double tolerance, int maxIterations) : NewtonRaphsonSolver<LegendreGrid1D>(actualClass, tolerance, maxIterations), actualClass(actualClass) {}

        double objectiveFunction(double x) {
            return actualClass.poly.evalD012(x)[0];
        }

        double objectiveDerivative(double x) {
            return actualClass.poly.evalD012(x)[1];
        }

    };

    LegendreGrid1D(std::size_t N) : order_(N), points_(N), weights_(N), poly(N) {

        RootSolver solver(*this, 1e-16, 1000);
        ChebyshevGrid1D<FloatingDataType> chebyGrid(N);
        for (int i = 0; i <= points_.size(); i++) {
            double xLower = chebyGrid.getPoints()[i] - 1;
            double xUpper = chebyGrid.getPoints()[i] + 1;
            points_[i] = solver.solve(xLower, xUpper);
        }

        // int p = order_;
        // int ph = floor((p + 1) / 2);

        // for (int i = 0; i < ph; i++) {
        //     double ii = i + 1.0;
        //     double x = cos((2.0*ii - 1.0)*M_PI / (2.0*p + 1.0));
        //     std::vector<double> L;
        //     for (int k = 0; k < 20; k++) {
        //         L = poly(x);
        //         double dx = -L[0] / L[1];
        //         x = x + dx;
        //         if (abs(x) < 1e-16) {
        //             break;
        //         }
        //     }
        //     points_[p + 2 - i] = x;
        //     weights_[p + 2 - i] = 2.0 / ((1.0 - pow(x,2))*pow(L[1],2));
        // }

        // if (p + 1 != 2*ph) {
        //     double x = 0.0;
        //     std::vector<double> L = poly(x);
        //     points_[ph + 1] = x;
        //     weights_[ph + 1] = 2.0 / ((1.0 - pow(x,2))*pow(L[1],2));
        // }

        // for (int i = 0; i < ph; i++) {
        //     points_[i] = -points_[p + 2 - i];
        //     weights_[i] = weights_[p + 2 - i];
        // }

    }

    std::vector<FloatingDataType> getPoints() { return points_; }
    std::size_t getNPoints() { return points_.size(); }
    FloatingDataType getLowerBound() { return -1.0; }
    FloatingDataType getUpperBound() { return 1.0; }

private:

    std::size_t order_;
    std::vector<FloatingDataType> points_;
    std::vector<FloatingDataType> weights_;
    

};

} // NAMESPACE: HydroForest

#endif // GRID_1D_HPP_