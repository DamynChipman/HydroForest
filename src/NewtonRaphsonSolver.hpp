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
            std::cerr << "[HydroForest::NewtonRaphsonSolver::solve] Root must be bracketed by `xLowerBound` and `xUpperBound`" << std::endl;
            throw HydroForestException("Root must be bracketed by `xLowerBound` and `xUpperBound`");
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
        std::cerr << "[HydroForest::NewtonRaphsonSolver::solve] Maximum number of iterations reached in solve!" << std::endl;
        throw HydroForestException("Maximum number of iterations reached in solve!");
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
        std::cerr << "[HydroForest::NewtonRaphsonSolver::solve] Maximum number of iterations reached in solve!" << std::endl;
        throw HydroForestException("Maximum number of iterations reached in solve!");
    }

};

};

#endif // NEWTON_RAPHSON_SOLVER_HPP_