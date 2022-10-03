#ifndef OPTIONS_HPP_
#define OPTIONS_HPP_

#include <iostream>
#include <string>
#include <map>
#include <variant>

namespace HydroForest {

/**
 * @brief Handles any user provided and simulation options
 * 
 */
class Options {

public:
    
    /**
     * @brief Set of possible option types
     * 
     * [std::string, bool, int, double]
     * 
     */
    using OptionTypes = std::variant<std::string, bool, int, double>;

    /**
     * @brief Construct a new Options object
     * 
     */
    Options() {}

    /**
     * @brief Construct a new Options object from an already built map
     * 
     * @param map 
     */
    Options(std::map<std::string, OptionTypes> map) : optionsMap_(map) {}

    /**
     * @brief Index into options map to get or set option
     * 
     * @param key Name of option
     * @return OptionTypes& Value of option
     */
    OptionTypes& operator[](std::string const key) {
        return optionsMap_[key];
    }

    /**
     * @brief Set the options from a given file
     * 
     * NOTE: Not implemented!
     * TODO: Implement!
     * 
     * @param filename 
     */
    void setFromFile(std::string filename) {

    }

    /**
     * @brief Print options to `os`
     * 
     * @param os Where to print to
     * @param options Options instance
     * @return std::ostream& ostream reference
     */
    friend std::ostream& operator<<(std::ostream& os, const Options& options) {
        for (const auto& [key, value] : options.optionsMap_) {
            os << key << " : ";
            std::visit([&] (auto&& v)
                {os << v << std::endl; }, value);
        }
        return os;
    }

private:

    /**
     * @brief Stores the options in a std::map
     * 
     */
    std::map<std::string, OptionTypes> optionsMap_;

};

} // NAMESPACE: HydroForest

#endif // OPTIONS_HPP_