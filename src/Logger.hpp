#ifndef LOGGER_HPP_
#define LOGGER_HPP_

#include <stdio.h>
#include <stdarg.h>
#include <string>
#include <mpi.h>

namespace HydroForest {

/**
 * @brief Struct to control output of messages to the console or a file
 * 
 */
struct Logger {

    /**
     * @brief Construct a new Logger object
     * 
     */
    Logger() {}

    /**
     * @brief Destroy the Logger object
     * 
     */
    ~Logger() {}

    /**
     * @brief Prints a message to the console with optional formatting for `printf`
     * 
     * @tparam Args 
     * @param message The message to log
     * @param args Additional arguments to pass to `printf` for formatting
     */
    template<class... Args>
    void log(std::string message, Args... args) {
        int myRank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        std::string toPrint = "[HydroForest " + std::to_string(myRank) + "] " + message + "\n";
        printf(toPrint.c_str(), args...);
    }

};

} // NAMESPACE: HydroForest

#endif // LOGGER_HPP_