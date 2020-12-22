#include <iostream>
#include "Solver.h"
int main() {
    std::cout << "Hello, World!" << std::endl;
    DGTCM::Solver solver(0.01, 0.1, 30.0);
    solver.Solve();
    return 0;
}
