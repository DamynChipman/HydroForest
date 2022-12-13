#include "P4est.hpp"

namespace HydroForest {

namespace p4est {

int p4est_to_HydroForest_face_index(int p4est_face_index) {
    switch (p4est_face_index) {
        case 0: return 3;
        case 1: return 1;
        case 2: return 0;
        case 3: return 2;
    }
}

int HydroForest_to_p4est_face_index(int HydroForest_face_index) {
    switch (HydroForest_face_index) {
        case 0: return 2;
        case 1: return 1;
        case 2: return 3;
        case 3: return 0;
    }
}

p4est_connectivity_t* p4est_connectivity_new_rectangular_domain(double x_lower, double x_upper, double y_lower, double y_upper) {
  
    const p4est_topidx_t num_vertices = 4;
    const p4est_topidx_t num_trees = 1;
    const p4est_topidx_t num_corners = 1;
    const double        vertices[4 * 3] = {
        x_lower, y_lower, 0,
        x_upper, y_lower, 0,
        x_lower, y_upper, 0,
        x_upper, y_upper, 0,
    };
    const p4est_topidx_t tree_to_vertex[1 * 4] = {
        0, 1, 2, 3
    };
    const p4est_topidx_t tree_to_tree[1 * 4] = {
        0, 0, 0, 0
    };
    const int8_t        tree_to_face[1 * 4] = {
        1, 0, 3, 2
        // 0, 1, 2, 3
    };
    const p4est_topidx_t tree_to_corner[1 * 4] = {
        0, 0, 0, 0
    };
    const p4est_topidx_t ctt_offset[1 + 1] = {
        0, 4
    };
    const p4est_topidx_t corner_to_tree[4] = {
        0, 0, 0, 0
    };
    const int8_t corner_to_corner[4] = {
        0, 1, 2, 3
    };

    return p4est_connectivity_new_copy (num_vertices, num_trees, num_corners,
                                        vertices, tree_to_vertex,
                                        tree_to_tree, tree_to_face,
                                        tree_to_corner, ctt_offset,
                                        corner_to_tree, corner_to_corner);
}
    // return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
    //                                     vertices, tree_to_vertex,
    //                                     tree_to_tree, tree_to_face,
    //                                     NULL, &ctt_offset[0],
    //                                     NULL, NULL);

} // NAMESPACE : p4est

} // NAMESPACE : HydroForest