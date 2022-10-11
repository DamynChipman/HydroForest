#ifndef TIME_INTEGRATION_HPP_
#define TIME_INTEGRATION_HPP_

#include <functional>

#include "Vector.hpp"
#include "Matrix.hpp"

namespace HydroForest {

class TimeIntegration {

public:

    virtual Vector<double> update(double t, double dt, Vector<double> q_n, Matrix<double> Rq_n) = 0;
    virtual Vector<double> update(double t, double dt, Vector<double> q_n, std::function<Vector<double>(Vector<double>)> Rq_n) = 0;

};

class RungeKutta3 : public TimeIntegration {

    RungeKutta3() {}

};

class ExplicitEuler : public TimeIntegration {

    ExplicitEuler() {}

};

} // NAMESPACE : HydroForest

#endif // TIME_INTEGRATION_HPP_