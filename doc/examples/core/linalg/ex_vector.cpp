
#include <iostream>
#include "simol/core/linalg/Vector.hpp"

int main()
{
    simol::Vector<double> v(3);
    v(0) = v(1) = v(2) = 1;
    std::cout << v << std::endl;
}
