#include "utilities/heterogeneous_containers.hpp"

#include <algorithm>
#include <iostream>
#include <string>

int main() {
    auto vec_writer = [](auto vec) {
        std::for_each(vec.cbegin(), vec.cend(), [](auto val) { std::cout << " " << val; });
    };

    auto map_writer = [](auto map) {
        std::for_each(map.cbegin(), map.cend(), [](auto val) { std::cout << " " << val.second; });
    };

    {  // testing vector functionality
        Utilities::HeterogeneousVector<double, int, std::string> vec;

        vec.template emplace_back<double>(1.);
        vec.template emplace_back<double>(2.);

        vec.template emplace_back<int>(1);
        vec.template emplace_back<int>(2);

        vec.template emplace_back<std::string>("Foo");

        if (vec.size() != 5) {
            return 1;
        }

        if (vec.at<double>(0) != 1.) {
            return 1;
        }

        if (vec.at<int>(0) != 1) {
            return 1;
        }

        if (vec.at<std::string>(0) != "Foo") {
            return 1;
        }

        std::cout << "We found " << vec.size() << "/5 elements in vec\n";
        std::cout << "They are: \n";

        Utilities::for_each_in_tuple(vec.data, vec_writer);

        std::cout << "\n";
    }

    {  // testing map functionality
        Utilities::HeterogeneousMap<double, int, std::string> map;

        map.template emplace<double>(0, 1.);
        map.template emplace<double>(10, 2.);

        map.template emplace<int>(5, 1);
        map.template emplace<int>(15, 2);

        map.template emplace<std::string>(20, "Foo");

        if (map.size() != 5) {
            return 1;
        }

        if (map.at<double>(0) != 1.) {
            return 1;
        }

        if (map.at<int>(5) != 1) {
            return 1;
        }

        if (map.at<std::string>(20) != "Foo") {
            return 1;
        }

        std::cout << "We found " << map.size() << "/5 elements in map\n";
        std::cout << "They are: \n";

        Utilities::for_each_in_tuple(map.data, map_writer);

        std::cout << "\n";
    }

    return 0;
}