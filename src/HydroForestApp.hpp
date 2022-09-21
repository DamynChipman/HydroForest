#ifndef HYDRO_FOREST_APP_HPP_
#define HYDRO_FOREST_APP_HPP_

#include <p4est.h>
#include <petsc.h>
#include "GenericSingleton.hpp"
#include "Logger.hpp"
#include "Options.hpp"

namespace HydroForest {

class HydroForestException : public std::exception {

public:

    explicit HydroForestException(const char* message) : message_(message) {}
    explicit HydroForestException(const std::string message) : message_(message) {}
    virtual ~HydroForestException() noexcept {}

    virtual const char* what() const noexcept {
        std::string w = "[HydroForest Exception] " + message_;
        return w.c_str();
    }

private:

    std::string message_;

};

class HydroForestApp : public GenericSingleton<HydroForestApp> {

private:

    int* argc_;
    char*** argv_;
    Logger logger_;
    Options options_;

public:

    HydroForestApp() : argc_(nullptr), argv_(nullptr) {}
    
    HydroForestApp(int* argc, char*** argv) : argc_(argc), argv_(argv) {
        MPI_Init(argc_, argv_);
        PetscInitialize(argc_, argv_, NULL, NULL);
        
        int myRank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        if (myRank == 0) {
            std::cout << "[HydroForest] Welcome to HydroForest!" << std::endl;
        }
        this->actualClassPointer_ = this;

    }

    ~HydroForestApp() {
        int myRank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        if (myRank == 0) {
            std::cout << "[HydroForest] End of app life cycle, finalizing..." << std::endl;
        }
        PetscFinalize();
        MPI_Finalize();
    }

    int* getArgc() const { return argc_; }
    char*** getArgv() const { return argv_; }
    Logger& getLogger() { return logger_; }
    Options& getOptions() { return options_; }

    template<class... Args>
    void log(std::string message, Args... args) {
        logger_.log(message, args...);
    }

};


} // NAMESPACE: HydroForest

#endif // HYDRO_FOREST_APP_HPP_