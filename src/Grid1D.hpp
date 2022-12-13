#ifndef GRID_1D_HPP_
#define GRID_1D_HPP_

#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <iostream>

#include "Vector.hpp"
#include "Polynomial.hpp"
#include "VTKBuilder.hpp"
#include "NewtonRaphsonSolver.hpp"

namespace HydroForest {

template<typename FloatingDataType>
class Grid1DBase {

public:

    virtual Vector<FloatingDataType>& getPoints() = 0;
    virtual Vector<FloatingDataType> getWeights() = 0;
    virtual std::size_t getNPoints() = 0;
    virtual FloatingDataType getLowerBound() = 0;
    virtual FloatingDataType getUpperBound() = 0;
    virtual FloatingDataType operator[](std::size_t index) = 0;

    friend std::ostream& operator<<(std::ostream& os, Grid1DBase& grid) {
        // os << "--- Grid ---" << std::endl;
        os << "Lower Bound = " << grid.getLowerBound() << std::endl;
        os << "Upper Bound = " << grid.getUpperBound() << std::endl;
        os << "# of Points = " << grid.getNPoints() << std::endl;
        for (auto i = 0; i < grid.getPoints().size(); i++) {
            os << std::setprecision(4) << std::setw(8) << grid.operator[](i);
            if (i % 10 == 9) os << std::endl;
        }
        os << std::endl;
        return os;
    }

};

template<typename FloatingDataType>
class UniformGrid1D : public Grid1DBase<FloatingDataType>, public RectilinearGridNodeBase, public DataArrayNodeBase {

public:

    UniformGrid1D(FloatingDataType lowerBound, FloatingDataType upperBound, std::size_t nPoints) : lowerBound_(lowerBound), upperBound_(upperBound), nPoints_(nPoints), points_(nPoints) {
        gridSpacing_ = (upperBound_ - lowerBound_) / (nPoints_ - 1);
        for (auto i = 0; i < nPoints_; i++) {
            points_[i] = lowerBound_ + i * gridSpacing_;
        }
    }

    Vector<FloatingDataType>& getPoints() { return points_; }
    Vector<FloatingDataType> getWeights() { return Vector<FloatingDataType>(points_.size(), 1.0); }
    std::size_t getNPoints() { return points_.size(); }
    FloatingDataType getLowerBound() { return lowerBound_; }
    FloatingDataType getUpperBound() { return upperBound_; }
    FloatingDataType operator[](std::size_t index) { return points_[index]; }

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
            // for (auto& p : points_) pointsAsString += std::to_string(p) + " ";
            return pointsAsString;
        }
    }

private:

    FloatingDataType lowerBound_;
    FloatingDataType upperBound_;
    FloatingDataType gridSpacing_;
    std::size_t nPoints_;
    Vector<FloatingDataType> points_;

};

template <typename FloatingDataType>
class ChebyshevGrid1D : public Grid1DBase<FloatingDataType> {

public:

    ChebyshevGrid1D(std::size_t order) : order_(order), points_(order+1), weights_(order+1) {

        for (int i = 0; i <= order; i++) {
            int index = order - i;
            points_[index] = cos((2.0*i + 1.0)/(2.0*order + 2.0) * M_PI);
            weights_[i] = (M_PI) / (order + 1.0);
        }

        // std::sort(points_.data().begin(), points_.data().end(), std::greater<FloatingDataType>{});
        // std::sort(weights_.data().begin(), weights_.data().end(), std::greater<FloatingDataType>{});

    }

    double weightFunction(double x) {
        return (1.0) / (sqrt(1.0 - pow(x, 2)));
    }

    Vector<FloatingDataType>& getPoints() { return points_; }
    Vector<FloatingDataType> getWeights() { return weights_; }
    std::size_t getNPoints() { return points_.size(); }
    FloatingDataType getLowerBound() { return -1.0; }
    FloatingDataType getUpperBound() { return 1.0; }
    FloatingDataType operator[](std::size_t index) { return points_[index]; }

private:

    std::size_t order_;
    Vector<FloatingDataType> points_;
    Vector<FloatingDataType> weights_;

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
            return this->actualClass.poly.derivatives012(x)[1];
        }

    };

    LegendreGrid1D(std::size_t order) : order_(order), points_(order+1), weights_(order+1), poly(order+1) {

        RootSolver solver(*this);
        solver.tolerance = 1e-15;
        ChebyshevGrid1D<FloatingDataType> chebyGrid(order);
        for (int i = 0; i < points_.size(); i++) {
            points_[i] = solver.solve(chebyGrid.getPoints()[i]);
            
            double dPhi = poly.derivatives012(points_[i])[1];
            weights_[i] = (2.0) / ((1.0 - pow(points_[i],2)) * pow(dPhi, 2));
        }

    }

    Vector<FloatingDataType>& getPoints() { return points_; }
    Vector<FloatingDataType> getWeights() { return weights_; }
    std::size_t getNPoints() { return points_.size(); }
    FloatingDataType getLowerBound() { return -1.0; }
    FloatingDataType getUpperBound() { return 1.0; }
    FloatingDataType operator[](std::size_t index) { return points_[index]; }

private:

    std::size_t order_;
    Vector<FloatingDataType> points_;
    Vector<FloatingDataType> weights_;

};

template <typename FloatingDataType>
class LobattoGrid1D : public Grid1DBase<FloatingDataType>, public DataArrayNodeBase {

public:

    LegendrePolynomial poly;

    struct RootSolver : public NewtonRaphsonSolver<LobattoGrid1D> {

        RootSolver(LobattoGrid1D& actualClass) : NewtonRaphsonSolver<LobattoGrid1D>(actualClass) {}

        double objectiveFunction(double x) {
            Vector<double> phiLegendreValues = this->actualClass.poly.derivatives012(x);
            return (1.0 - pow(x,2))*phiLegendreValues[1];
        }

        double objectiveDerivative(double x) {
            Vector<double> phiLegendreValues = this->actualClass.poly.derivatives012(x);
            return -2.0*x*phiLegendreValues[1] + (1.0 - pow(x,2))*phiLegendreValues[2];
        }

    };

    LobattoGrid1D(std::size_t order) : order_(order), points_(order+1), weights_(order+1), poly(order) {

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

    std::size_t order() { return order_; }
    Vector<FloatingDataType>& getPoints() { return points_; }
    Vector<FloatingDataType> getWeights() { return weights_; }
    std::size_t getNPoints() { return points_.size(); }
    FloatingDataType getLowerBound() { return -1.0; }
    FloatingDataType getUpperBound() { return 1.0; }
    FloatingDataType operator[](std::size_t index) { return points_[index]; }

    std::string getType() {
        return "Float32";
    }

    std::string getName() {
        return "LobattoGrid";
    }

    std::string getNumberOfComponents() {
        return "1";
    }

    std::string getFormat() {
        return "ascii";
    }

    std::string getRangeMin() {
        return "-1.0";
    }

    std::string getRangeMax() {
        return "1.0";
    }

    std::string getData() {
        std::string pointsAsString = "";
        for (auto i = 0; i < points_.size(); i++) {
            pointsAsString+= std::to_string(points_[i]) + " ";
        }
        // for (auto& p : points_) pointsAsString += std::to_string(p) + " ";
        return pointsAsString;
    }


private:

    std::size_t order_;
    Vector<FloatingDataType> points_;
    Vector<FloatingDataType> weights_;

};

} // NAMESPACE: HydroForest

#endif // GRID_1D_HPP_