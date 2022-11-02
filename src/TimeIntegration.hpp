#ifndef TIME_INTEGRATION_HPP_
#define TIME_INTEGRATION_HPP_

#include <functional>

#include "Vector.hpp"
#include "Matrix.hpp"

namespace HydroForest {

template<typename FloatingDataType>
class TimeIntegration {

public:

    virtual Vector<FloatingDataType> update(FloatingDataType t, FloatingDataType dt, Vector<FloatingDataType> q_n, Matrix<FloatingDataType> Rq_n) = 0;
    virtual Vector<FloatingDataType> update(FloatingDataType t, FloatingDataType dt, Vector<FloatingDataType> q_n, Vector<FloatingDataType> Rq_n) = 0;
    virtual Vector<FloatingDataType> update(FloatingDataType t, FloatingDataType dt, Vector<FloatingDataType> q_n, std::function<Vector<FloatingDataType>(Vector<FloatingDataType>)> Rq_n) = 0;
    virtual FloatingDataType getCFLStabilityNumber() = 0;
    virtual FloatingDataType getMaxTimeStep(FloatingDataType characteristicSpeed, FloatingDataType characteristicLength) = 0;

};

template<typename FloatingDataType>
class RungeKutta3 : public TimeIntegration<FloatingDataType> {

public:

    RungeKutta3() {}

    Vector<FloatingDataType> update(FloatingDataType t, FloatingDataType dt, Vector<FloatingDataType> q_n, Matrix<FloatingDataType> Rq_n) {

        Vector<FloatingDataType> q_temp;
        q_temp = Rq_n*q_n;
        Vector<FloatingDataType> q_1 = q_n + dt*q_temp;
        q_temp = Rq_n*q_1;
        Vector<FloatingDataType> q_2 = (3.0/4.0)*q_n + (1.0/4.0)*q_1 + (1.0/4.0)*dt*q_temp;
        q_temp = Rq_n*q_2;
        Vector<FloatingDataType> q_np1 = (1.0/3.0)*q_n + (2.0/3.0)*q_2 + (2.0/3.0)*dt*q_temp;
        return q_np1;

    }

    Vector<FloatingDataType> update(FloatingDataType t, FloatingDataType dt, Vector<FloatingDataType> q_n, std::function<Vector<FloatingDataType>(Vector<FloatingDataType>)> Rq_n) {

        Vector<FloatingDataType> q_1 = q_n + dt*Rq_n(q_n);
        Vector<FloatingDataType> q_2 = (3.0/4.0)*q_n + (1.0/4.0)*q_1 + (1.0/4.0)*dt*Rq_n(q_1);
        Vector<FloatingDataType> q_np1 = (1.0/3.0)*q_n + (2.0/3.0)*q_2 + (2.0/3.0)*dt*Rq_n(q_2);
        return q_np1;

    }

    Vector<FloatingDataType> update(FloatingDataType t, FloatingDataType dt, Vector<FloatingDataType> q_n, Vector<FloatingDataType> Rq_n) {

        Vector<FloatingDataType> q_1 = q_n + dt*Rq_n;
        Vector<FloatingDataType> q_2 = (3.0/4.0)*q_n + (1.0/4.0)*q_1 + (1.0/4.0)*dt*Rq_n;
        Vector<FloatingDataType> q_np1 = (1.0/3.0)*q_n + (2.0/3.0)*q_2 + (2.0/3.0)*dt*Rq_n;
        return q_np1;

    }

    FloatingDataType getCFLStabilityNumber() { return cfl_; }
    FloatingDataType getMaxTimeStep(FloatingDataType characteristicSpeed, FloatingDataType characteristicLength) {
        return (characteristicLength*cfl_) / characteristicSpeed;
    }

protected:

    FloatingDataType cfl_ = 1.0/3.0;

};

template<typename FloatingDataType>
class ExplicitEuler : public TimeIntegration<FloatingDataType> {

    ExplicitEuler() {}

};

} // NAMESPACE : HydroForest

#endif // TIME_INTEGRATION_HPP_