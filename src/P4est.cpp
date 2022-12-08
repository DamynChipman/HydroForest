#include "P4est.hpp"

namespace HydroForest {

namespace p4est {

p4est_connectivity_t* p4est_connectivity_new_rectangular_domain(double x_lower, double x_upper, double y_lower, double y_upper) {
  
    const p4est_topidx_t num_vertices = 4;
    const p4est_topidx_t num_trees = 1;
    const p4est_topidx_t num_ctt = 0;
    const double        vertices[4 * 3] = {
        x_lower, y_lower, 0,
        x_upper, y_lower, 0,
        x_lower, y_upper, 0,
        x_upper, y_upper, 0,
    };
    const p4est_topidx_t tree_to_vertex[1 * 4] = {
        0, 1, 2, 3,
    };
    const p4est_topidx_t tree_to_tree[1 * 4] = {
        0, 0, 0, 0,
    };
    const int8_t        tree_to_face[1 * 4] = {
        0, 1, 2, 3,
    };

    return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                        vertices, tree_to_vertex,
                                        tree_to_tree, tree_to_face,
                                        NULL, &num_ctt, NULL, NULL);

}

} // NAMESPACE : p4est

} // NAMESPACE : HydroForest