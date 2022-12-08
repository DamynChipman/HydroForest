#ifndef TIMER_HPP_
#define TIMER_HPP_

#include <chrono>

namespace HydroForest {

class Timer {

public:

    using Clock = std::chrono::steady_clock;

    Timer() {}

    void start() {
        startTime_ = Clock::now();
    }

    double stop() {
        endTime_ = Clock::now();
        std::chrono::duration<double> diff = endTime_ - startTime_;
        accumulatedTime_ += diff.count();
        return diff.count();
    }

    double time() { return accumulatedTime_; }

    void restart() {
        accumulatedTime_ = 0;
    }

private:

    double accumulatedTime_ = 0;
    std::chrono::time_point<Clock> startTime_;
    std::chrono::time_point<Clock> endTime_;

};

} // NAMESPACE : HydroForest

#endif // TIMER_HPP_