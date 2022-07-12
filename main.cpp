#include <iostream>

#include "kylin/matrix.hpp"
using namespace kylin;

int main() {
    DynamicMatrix<int> A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    std::cout << A << std::endl;
    StaticMatrix<int, 3, 3> B = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    std::cout << B << std::endl;

    std::cout << (A + B) << std::endl;
    std::cout << (A - B) << std::endl;
    std::cout << (A * B) << std::endl;
    std::cout << A.transpose() << std::endl;

    std::cout << diag<int>({1, 2, 3}) << std::endl;
    std::cout << eye<int>(3, 1) << std::endl;
    return 0;
}
