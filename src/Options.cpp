#include "Options.hpp"

namespace HydroForest {

Options::OptionTypes Options::operator[](std::string const key) {
    return optionsMap_[key];
}

void Options::setOption(std::string const key, OptionTypes const value) {
    optionsMap_[key] = value;
}

void Options::setFromFile(std::string filename) {

}

std::ostream& operator<<(std::ostream& os, const Options& options) {
    for (const auto& [key, value] : options.optionsMap_) {
        os << key << " : ";
        std::visit([&] (auto&& v)
            {os << v << std::endl; }, value);
    }
    return os;
}


} // NAMESPACE: HydroForest