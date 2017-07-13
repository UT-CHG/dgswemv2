#include "utilities/heterogeneous_containers.hpp"

#include <algorithm>
#include <iostream>
#include <string>

int main() {

  auto writer = [](auto vec) {
    std::for_each(vec.cbegin(), vec.cend(),[](auto val) {
        std::cout << " " << val;
      });
  };

  {// testing generic functionality
    Utilities::HeterogeneousVector<double,int,std::string> vec;

    vec.template emplace_back<double>(1.);
    vec.template emplace_back<double>(2.);
    vec.template emplace_back<double>(2.);
    vec.template emplace_back<double>(2.);

    vec.template emplace_back<std::string>("Foo");

    std::cout << "We found " << vec.size() << "/5 elements in vec\n";
    std::cout << "They are: \n";

    if ( vec.size() != 5 ) {
      return 1;
    }

    Utilities::for_each_in_tuple(vec.data, writer);

    std::cout << "\n";
  }

  {//testing case when we only have one type
    Utilities::HeterogeneousVector<double> n_vec;

    n_vec.template emplace_back<double>(1.);
    n_vec.template emplace_back<double>(1.);
    n_vec.template emplace_back<double>(2.);
    n_vec.template emplace_back<double>(2.);
    n_vec.template emplace_back<double>(2.);

    std::cout << "We found " << n_vec.size() << "/5 elements in n_vec\n";
    std::cout << "They are: \n";


    Utilities::for_each_in_tuple(n_vec.data, writer);

    std::cout << "\n";

    if ( n_vec.size() != 5 ) {
      return 1;
    }
  }

  return 0;
}
