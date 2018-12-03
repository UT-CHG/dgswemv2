#include <cassert>
#include <iostream>
#include <tuple>
#include <vector>

#include "utilities/is_soa.hpp"


struct Accessor {
    int x;
    int y;
};

struct SoA {
public:
    using AccessorType = Accessor;

    Accessor at(unsigned int index) {
        assert(index < x.size() && index < y.size());
    }

private:
    std::vector<int> x;
    std::vector<int> y;
};


using NotSoA = std::tuple<>;

//has AccessorType but no at operation
struct MalformedSoA {
    using AccessorType = Accessor;
};

int main() {
    bool error_found{false};

    { // check SoA
        bool is_SoA_val = Utilities::is_SoA<SoA>::value;
        std::cout << std::boolalpha
                  << "Utilities::is_SoA<SoA>::value = " << is_SoA_val
                  << " (should be true)\n\n";

        //if is_SoA_val is false return error
        error_found |= !is_SoA_val;
    }

    { // check NotSoA
        bool is_SoA_val = Utilities::is_SoA<NotSoA>::value;
        std::cout << std::boolalpha
                  << "Utilities::is_SoA<NotSoA>::value = " << is_SoA_val
                  << " (should be false)\n\n";

        //if is_SoA_val is true return error
        error_found |= is_SoA_val;
    }

    { // vector has at() but should not be considered an SoA
        bool is_SoA_val = Utilities::is_SoA<std::vector<int>>::value;
        std::cout << std::boolalpha
                  << "Utilities::is_SoA<std::vector<int>>::value = " <<is_SoA_val
                  << " (should be false)\n\n";

        //if is_SoA_val is true return error
        error_found |= is_SoA_val;
    }

    { // malformed accessor (no at() member function
        bool is_SoA_val = Utilities::is_SoA<MalformedSoA>::value;
        std::cout << std::boolalpha
                  << "Utilities::is_SoA<MalformedSoA>::value = " <<is_SoA_val
                  << " (should be false)\n\n";

        //if is_SoA_val is true return error
        error_found |= is_SoA_val;
    }

    return error_found;
}