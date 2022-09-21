#ifndef NODAL_POINTS_HPP_
#define NODAL_POINTS_HPP_

#include <vector>
#include <math.h>

namespace HydroForest {

class NodalPointsBase {

public:



};

class ChebyshevPointsBasis {

public:

    ChebyshevPointsBasis(int order) : order_(order), points_(order) {

        for (int i = 0; i < order_; i++) {
            points_[i] = cos((2.0*i + 1.0)/(2.0*order_ + 2.0) * M_PI);
        }

    }

    double operator()(double x) {

    }

private:

    int order_;
    std::vector<double> points_;

};

} // NAMESPACE : HydroForest

#endif // NODAL_POINTS_HPP_