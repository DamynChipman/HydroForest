#ifndef OPTIONS_HPP_
#define OPTIONS_HPP_

#include <iostream>
#include <string>
#include <map>
#include <variant>

namespace HydroForest {

class Options {

public:
    
    using OptionTypes = std::variant<std::string, bool, int, double>;

    Options() {}
    Options(std::map<std::string, OptionTypes> map) : optionsMap_(map) {}

    OptionTypes operator[](std::string const key);
    void setOption(std::string const key, OptionTypes const value);
    void setFromFile(std::string filename);
    friend std::ostream& operator<<(std::ostream& os, const Options& options);

private:

    std::map<std::string, OptionTypes> optionsMap_;

};

} // NAMESPACE: HydroForest

#endif // OPTIONS_HPP_