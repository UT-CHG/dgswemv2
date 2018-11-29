#include "utilities/at_each.hpp"

#include <iostream>
#include <vector>

struct Accessor {
    Accessor()=delete;
    int& x;
    int& y;
};

struct AoS {
public:
    using AccessorType = Accessor;

    AoS()=default;
    AoS(size_t sz, size_t offset) {
        x.reserve(sz);
        y.reserve(sz);

        for (uint i = 0; i < sz; ++i ) {
            x.push_back(i + offset);
            y.push_back(i + sz + offset);
        }
    };

    AccessorType at(uint index) {
        return AccessorType{x.at(index), y.at(index)};
    }

private:
    std::vector<int> x;
    std::vector<int> y;
};


int main() {
    bool error_found = false;

    { //check that the AoS is working properly
        std::vector<AoS> v_aos;
        v_aos.reserve(3);
        v_aos.emplace_back(10,0);
        v_aos.emplace_back(10,1);
        v_aos.emplace_back(10,2);

        std::vector<Accessor> v_acsr = Utilities::at_each(v_aos,2);

        std::cout << "vector of AoS\n"
                  << " Computed    |   Expected\n"
                  << "---------------------------------------------\n"
                  << ' ' << v_acsr[0].x << "  |  " << 2 << '\n'
                  << ' ' << v_acsr[1].x << "  |  " << 3 << '\n'
                  << ' ' << v_acsr[2].x << "  |  " << 4 << '\n'
                  << v_acsr[0].y << "  |  " << 12<< '\n'
                  << v_acsr[1].y << "  |  " << 13<< '\n'
                  << v_acsr[2].y << "  |  " << 14<< '\n';

        error_found |= !( (v_acsr[0].x == 2) && (v_acsr[0].y == 12 ) &&
                          (v_acsr[1].x == 3) && (v_acsr[1].y == 13 ) &&
                          (v_acsr[2].x == 4) && (v_acsr[2].y == 14 ) );
    }

    { //check whether it works with a vector of vectors
        std::vector<std::vector<int>> v_vec{ std::vector<int>{0,4,5},
                                             std::vector<int>{1,2,3},
                                             std::vector<int>{9,8,7} };

        std::vector<std::reference_wrapper<int>> v_ref = Utilities::at_each(v_vec, 2);

        error_found |= !( (v_ref[0] == 5) &&
                          (v_ref[1] == 3) &&
                          (v_ref[2] == 7) );
    }

    { //check whether it works with an array of AoS
        std::array<AoS,3> ar_aos { AoS(10,0), AoS(10,1), AoS(10,2) };

        std::array<Accessor,3> ar_acsr = Utilities::at_each(ar_aos,2);

        std::cout << "array of AoS\n"
                  << " Computed    |   Expected\n"
                  << "---------------------------------------------\n"
                  << ' ' << ar_acsr[0].x << "  |  " << 2 << '\n'
                  << ' ' << ar_acsr[1].x << "  |  " << 3 << '\n'
                  << ' ' << ar_acsr[2].x << "  |  " << 4 << '\n'
                  << ar_acsr[0].y << "  |  " << 12<< '\n'
                  << ar_acsr[1].y << "  |  " << 13<< '\n'
                  << ar_acsr[2].y << "  |  " << 14<< '\n';

        error_found |= !( (ar_acsr[0].x == 2) && (ar_acsr[0].y == 12 ) &&
                          (ar_acsr[1].x == 3) && (ar_acsr[1].y == 13 ) &&
                          (ar_acsr[2].x == 4) && (ar_acsr[2].y == 14 ) );

    }

    return error_found;
}