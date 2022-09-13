#ifndef LOGGER_HPP_
#define LOGGER_HPP_

#include <stdio.h>
#include <stdarg.h>
#include <string>
#include <mpi.h>

namespace HydroForest {

struct Logger {

    Logger() {}
    ~Logger() {}

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