#ifndef SPACE_VECTOR_HPP_
#define SPACE_VECTOR_HPP_

#include "Vector.hpp"

namespace HydroForest {

template<typename FloatingDataType>
struct SpaceVector2D : public Vector<FloatingDataType> {

    SpaceVector2D(FloatingDataType x, FloatingDataType y) :
        Vector<FloatingDataType>({x,y})
            {}

    SpaceVector2D(Vector<FloatingDataType> v) :
        Vector<FloatingDataType>({v[0], v[1]})
            {}

    FloatingDataType x() { return this->operator[](0); }
    FloatingDataType y() { return this->operator[](1); }

};

template<typename FloatingDataType>
struct SpaceVector3D : public Vector<FloatingDataType> {

    SpaceVector3D(FloatingDataType x, FloatingDataType y, FloatingDataType z) :
        Vector<FloatingDataType>({x,y,z})
            {}

    SpaceVector3D(SpaceVector2D<FloatingDataType> v) :
        Vector<FloatingDataType>({v[0], v[1], 0})
            {}

    FloatingDataType x() { return this->operator[](0); }
    FloatingDataType y() { return this->operator[](1); }
    FloatingDataType z() { return this->operator[](2); }

};

} // NAMESPACE : HydroForest

#endif // SPACE_VECTOR_HPP_