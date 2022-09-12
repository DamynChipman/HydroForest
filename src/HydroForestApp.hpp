#ifndef HYDRO_FOREST_APP_HPP_
#define HYDRO_FOREST_APP_HPP_

#include <p4est.h>
#include <petsc.h>
#include "GenericSingleton.hpp"

namespace hf {

class HydroForestApp : public GenericSingleton<HydroForestApp> {

private:

    int* argc_;
    char*** argv_;

public:
    
    HydroForestApp(int* argc, char*** argv) : argc_(argc), argv_(argv) {
        std::cout << "[HydroForest] Welcome to HydroForest!" << std::endl;
        std::cout << "[HydroForest] Initializing MPI and PETSc..." << std::endl;
        MPI_Init(argc_, argv_);
        PetscInitialize(argc_, argv_, NULL, NULL);
    }

    ~HydroForestApp() {
        std::cout << "[HydroForest] End of app life cycle, finalizing..." << std::endl;
        PetscFinalize();
        MPI_Finalize();
    }

    int* getArgc() const { return argc_; }
    char*** getArgv() const { return argv_; }

};


} // NAMESPACE: hf

#endif // HYDRO_FOREST_APP_HPP_