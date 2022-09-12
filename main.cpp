#include <iostream>
#include <mpi.h>
#include <p4est.h>
#include <p4est_connectivity.h>
#include <p4est_vtk.h>
#include <petsc.h>

typedef int (*p4est_refine_t) (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant);

int refine(p4est* p4est, p4est_topidx_t which_tree, p4est_quadrant_t* quadrant)
{
    if (quadrant->level < 4)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int main(int argc, char** argv) {
    std::cout << "Hello, p4est!" << std::endl;

    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, NULL, NULL);

    p4est_connectivity_t* conn = p4est_connectivity_new_disk2d();
    p4est_geometry_t* geo = p4est_geometry_new_disk2d(conn, 1.0, 2.0);
    // p4est_connectivity_t* refinedConn = p4est_connectivity_refine(conn, 5);
    // p4est_t* p4est = p4est_new(MPI_COMM_WORLD, refinedConn, 0, NULL, NULL);
    p4est_t* p4est = p4est_new(MPI_COMM_WORLD, conn, 0, NULL, NULL);

    // p4est_refine(p4est, 1, &refine, NULL);

    p4est_vtk_write_file(p4est, geo, "p4est");
    // p4est_vtk_write_file(p4est, NULL, "p4est");
    p4est_destroy(p4est);
    // p4est_connectivity_destroy(conn);

    PetscFinalize();
    MPI_Finalize();
    
    return EXIT_SUCCESS;
}