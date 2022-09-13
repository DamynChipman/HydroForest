#include <iostream>

#include <HydroForestApp.hpp>
#include <UniformGrid1D.hpp>
#include <Options.hpp>

int main(int argc, char** argv) {

    HydroForest::HydroForestApp app(&argc, &argv);

    HydroForest::Options options;

    options.setOption("foo", "foo_value");
    options.setOption("bar", 42);
    options.setOption("pi", 3.14);

    std::cout << options << std::endl;

    return EXIT_SUCCESS;
}