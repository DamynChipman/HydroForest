#ifndef NEWTON_RAPHSON_SOLVER_HPP_
#define NEWTON_RAPHSON_SOLVER_HPP_

#include <functional>
#include <exception>
#include <string>

#include "HydroForestApp.hpp"

namespace HydroForest {

template<class ClassName>
class NewtonRaphsonSolver {

public:

    ClassName& actualClass;
    double tolerance = 1e-12;
    int maxIterations = 1000;

    virtual double objectiveFunction(double x) = 0;
    virtual double objectiveDerivative(double x) = 0;

    NewtonRaphsonSolver(ClassName& actualClass) : actualClass(actualClass), tolerance(1e-16), maxIterations(1000) {}
    NewtonRaphsonSolver(ClassName& actualClass, double tolerance, int maxIterations) : actualClass(actualClass), tolerance(tolerance), maxIterations(maxIterations) {}

    double solve(double xLowerBound, double xUpperBound) {
        double xLow, xHigh;
        double fLow = objectiveFunction(xLowerBound);
        double fHigh = objectiveFunction(xUpperBound);
        if ((fLow > 0.0 && fHigh > 0.0) || (fLow < 0.0 && fHigh < 0.0)) {
            std::string errorMessage = "[HydroForest::NewtonRaphsonSolver::solve] Root must be bracketed by `xLowerBound` and `xUpperBound`:\n";
            errorMessage += "\txLowerBound = " + std::to_string(xLowerBound) + ", xUpperBound = " + std::to_string(xUpperBound) + "\n";
            errorMessage += "\tfLow = " + std::to_string(fLow) + ", fHigh = " + std::to_string(fHigh) + "\n";
            std::cerr << errorMessage << std::endl;
            throw std::invalid_argument(errorMessage);
        }

        if (fLow == 0.0) { return xLowerBound; }
        if (fHigh == 0.0) { return xUpperBound; }
        if (fLow < 0.0) {
            xLow = xLowerBound;
            xHigh = xUpperBound;
        }
        else {
            xHigh = xLowerBound;
            xLow = xUpperBound;
        }
        double rts = 0.5*(xLowerBound + xUpperBound);
        double dxOld = abs(xUpperBound - xLowerBound);
        double dx = dxOld;
        double f = objectiveFunction(rts);
        double df = objectiveDerivative(rts);
        for (int j = 0; j < maxIterations; j++) {
            if ((((rts - xHigh)*df - f)*((rts - xLow)*df - f) > 0.0) || (abs(2.0*f) > abs(dxOld*df))) {
                dxOld = dx;
                dx = 0.5*(xHigh - xLow);
                rts = xLow + dx;
                if (xLow == rts) return rts;
            }
            else {
                dxOld = dx;
                dx = f / df;
                double temp = rts;
                rts -= dx;
                if (temp == rts) return rts;
            }
            if (abs(dx) < tolerance) return rts;
            double f = objectiveFunction(rts);
            double df = objectiveDerivative(rts);
            if (f < 0.0) {
                xLow = rts;
            }
            else {
                xHigh = rts;
            }
        }
        std::string errorMessage = "[HydroForest::NewtonRaphsonSolver::solve] Maximum number of iterations reached in solve:\n";
        errorMessage += "\tcurrent value = " + std::to_string(rts) + "\n";
        errorMessage += "\tmaxIterations = " + std::to_string(maxIterations) + "\n";
        errorMessage += "\ttolerance = " + std::to_string(tolerance) + "\n";
        std::cerr << errorMessage << std::endl;
        throw std::out_of_range(errorMessage);
    }

    double solve(double xGuess) {
        double x_j = xGuess;
        double x_jp1;
        for (int j = 0; j < maxIterations; j++) {
            double f = objectiveFunction(x_j);
            double df = objectiveDerivative(x_j);
            x_jp1 = x_j - (f / df);
            if (abs(x_jp1 - x_j) < tolerance) {
                return x_jp1;
            }
            x_j = x_jp1;
        }
        std::string errorMessage = "[HydroForest::NewtonRaphsonSolver::solve] Maximum number of iterations reached in solve:\n";
        errorMessage += "\tcurrent value = " + std::to_string(x_j) + "\n";
        errorMessage += "\tmaxIterations = " + std::to_string(maxIterations) + "\n";
        errorMessage += "\ttolerance = " + std::to_string(tolerance) + "\n";
        std::cerr << errorMessage << std::endl;
        throw std::out_of_range(errorMessage);
    }

};

};

#endif // NEWTON_RAPHSON_SOLVER_HPP_