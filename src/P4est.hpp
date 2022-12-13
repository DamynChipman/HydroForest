#ifndef P4EST_HPP_
#define P4EST_HPP_

#include <p4est.h>
#include <p4est_extended.h>
#include <p4est_connectivity.h>
#include <p4est_iterate.h>
#include <p4est_geometry.h>
#include <p4est_vtk.h>

namespace HydroForest {

namespace p4est {

int p4est_to_HydroForest_face_index(int p4est_face_index);
int HydroForest_to_p4est_face_index(int HydroForest_face_index);
p4est_connectivity_t* p4est_connectivity_new_rectangular_domain(double x_lower, double x_upper, double y_lower, double y_upper);

} // NAMESPACE : p4est

} // NAMESPACE : HydroForest

#endif // P4EST_HPP_