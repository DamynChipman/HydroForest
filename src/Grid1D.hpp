#ifndef GRID_1D_HPP_
#define GRID_1D_HPP_

#include <cmath>
#include <vector>
#include <algorithm>

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

        for (int i = 0; i <= order; i++) {
            points_[i] = cos((2.0*i + 1.0)/(2.0*order + 2.0) * M_PI);
            weights_[i] = (M_PI) / (order + 1.0);
        }

        std::sort(points_.begin(), points_.end());
        std::sort(weights_.begin(), weights_.end());

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

        // RootSolver() : NewtonRaphsonSolver<LegendreGrid1D>() {}
        RootSolver(LegendreGrid1D& actualClass) : NewtonRaphsonSolver<LegendreGrid1D>(actualClass) {}

        double objectiveFunction(double x) {
            double f = this->actualClass.poly(x);
            return this->actualClass.poly(x);
        }

        double objectiveDerivative(double x) {
            return this->actualClass.poly.evalD012(x)[1];
        }

    };

    LegendreGrid1D(std::size_t order) : order_(order), points_(order+1), weights_(order+1), poly(order+1) {

        RootSolver solver(*this);
        solver.tolerance = 1e-15;
        ChebyshevGrid1D<FloatingDataType> chebyGrid(order);
        for (int i = 0; i < points_.size(); i++) {
            points_[i] = solver.solve(chebyGrid.getPoints()[i]);
            
            double dPhi = poly.evalD012(points_[i])[1];
            weights_[i] = (2.0) / ((1.0 - pow(points_[i],2)) * pow(dPhi, 2));
        }

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
class LobattoGrid1D : public Grid1DBase<FloatingDataType> {

public:

    LegendrePolynomial poly;

    struct RootSolver : public NewtonRaphsonSolver<LobattoGrid1D> {

        RootSolver(LobattoGrid1D& actualClass) : NewtonRaphsonSolver<LobattoGrid1D>(actualClass) {}

        double objectiveFunction(double x) {
            std::vector<double> phiLegendreValues = this->actualClass.poly.evalD012(x);
            return (1.0 - pow(x,2))*phiLegendreValues[1];
        }

        double objectiveDerivative(double x) {
            std::vector<double> phiLegendreValues = this->actualClass.poly.evalD012(x);
            return -2.0*x*phiLegendreValues[1] + (1.0 - pow(x,2))*phiLegendreValues[2];
        }

    };

    LobattoGrid1D(std::size_t order) : order_(order), points_(order+1, 0), weights_(order+1), poly(order) {

        RootSolver solver(*this);
        solver.tolerance = 1e-15;
        ChebyshevGrid1D<FloatingDataType> chebyGrid(order);

        points_[0] = -1.0;
        points_[points_.size()-1] = 1.0;

        if (points_.size() % 2) {
            // Odd number of points; end points are -1, 1; middle is 0; rest are found via solver
            points_[(points_.size()-1) / 2] = 0.0;
            for (int i = 1; i < points_.size()-1; i++) {
                if (i != (points_.size()-1)/2) {
                    points_[i] = solver.solve(chebyGrid.getPoints()[i]);
                }
            }
        }
        else {
            // Even number of points; end points are -1, 1; rest are found via solver
            for (int i = 1; i < points_.size()-1; i++) {
                points_[i] = solver.solve(chebyGrid.getPoints()[i]);
            }
        }

        for (int i = 0; i < weights_.size(); i++) {
            weights_[i] = (2.0) / (order*(order + 1)*pow(poly(points_[i]),2));
        }

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