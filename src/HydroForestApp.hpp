#ifndef HYDRO_FOREST_APP_HPP_
#define HYDRO_FOREST_APP_HPP_

#include <iostream>
#include <string>
#include <map>
#include <variant>
#include <cstdarg>
#include <p4est.h>
#include <petsc.h>
#include "GenericSingleton.hpp"
#include "Timer.hpp"
#include "Logger.hpp"
#include "Options.hpp"
#include "Memory.hpp"

namespace HydroForest {

/**
 * @brief Single instance of the HydroFoest app
 * 
 * The user will create an instance of the HydroForest within the `main` function. It is derived
 * from a generic singleton design, meaning there is only one instance. At initialization, any
 * dependencies are initialized and the command line arguments are parsed and any HydroForest
 * options are stored in the `Options` class.
 * 
 */
class HydroForestApp : public GenericSingleton<HydroForestApp> {

protected:

    /**
     * @brief Address of command line arguments `argc`
     * 
     */
    int* argc_;

    /**
     * @brief Address of command line arguments `argv`
     * 
     */
    char*** argv_;

    /**
     * @brief Logger instance for logging during app
     * 
     */
    Logger logger_;

    /**
     * @brief Options instance for storing and setting options
     * 
     */
    Options options_;

    std::map<std::string, Timer> timers_;

public:

    /**
     * @brief Construct a new Hydro Forest App object
     * 
     * Default constructor; sets pointers to `nullptr`. The other constructor should be used
     * instead.
     * 
     */
    HydroForestApp() : argc_(nullptr), argv_(nullptr) {}
    
    /**
     * @brief Construct a new Hydro Forest App object from the command line arguments
     * 
     * Also initializes any dependencies:
     *      MPI
     *      PETSc
     * 
     * @param argc 
     * @param argv 
     */
    HydroForestApp(int* argc, char*** argv) : argc_(argc), argv_(argv) {
        MPI_Init(argc_, argv_);
        PetscInitialize(argc_, argv_, NULL, NULL);
        
        int myRank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        if (myRank == 0) {
            std::cout << "[HydroForest] Welcome to HydroForest!" << std::endl;
        }
        this->actualClassPointer_ = this;

        // Set options from command line
        setCMLOptions_();

    }

    /**
     * @brief Destroy the Hydro Forest App object
     * 
     */
    ~HydroForestApp() {
        int myRank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        if (myRank == 0) {
            std::cout << "[HydroForest] End of app life cycle, finalizing..." << std::endl;
        }
        PetscFinalize();
        MPI_Finalize();
    }

    /**
     * @brief Get the command line arguments `argc`
     * 
     * @return int* 
     */
    int* getArgc() const { return argc_; }

    /**
     * @brief Get the command line arguments `argv`
     * 
     * @return char*** 
     */
    char*** getArgv() const { return argv_; }

    /**
     * @brief Get the Logger instance
     * 
     * @return Logger& 
     */
    Logger& getLogger() { return logger_; }

    /**
     * @brief Get the Options instance
     * 
     * @return Options& 
     */
    Options& getOptions() { return options_; }

    /**
     * @brief Get the timers
     * 
     * @return std::map<std::string, Timer> 
     */
    std::map<std::string, Timer> timers() { return timers_; }

    /**
     * @brief Wrapper for `Logger.log()`
     * 
     * @tparam Args 
     * @param message The message to log
     * @param args Additional arguments to pass to `printf` for formatting
     */
    template<class... Args>
    void log(std::string message, Args... args) {
        logger_.log(message, args...);
    }

private:

    void setCMLOptions_() {
        // Loop through argv

    }

};

} // NAMESPACE: HydroForest

#endif // HYDRO_FOREST_APP_HPP_