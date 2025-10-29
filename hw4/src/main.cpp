#include <iomanip>
#include <iostream>
#include "Matrix.h"
#include "SystemEquation.h"
#include "solver/SORSolver.h"
#include "solver/LUSolver.h"

void test_matrix()
{
    tinyLinAlg::Matrix A(3, 3);
    A(0, 0) = 3;
    A(0, 1) = 2;
    A(0, 2) = -4;
    A(1, 0) = 2;
    A(1, 1) = 3;
    A(1, 2) = 3;
    A(2, 0) = 5;
    A(2, 1) = -3;
    A(2, 2) = 1;

    tinyLinAlg::Matrix B(3, 1);
    B(0, 0) = 3;
    B(1, 0) = 15;
    B(2, 0) = 14;

    A.print();
    B.print();
}

void test_system_equation()
{
    using namespace tinyLinAlg;

    // 1. 10×10 系数矩阵
    auto A = std::make_unique<Matrix>(10, 10);
    const double a[100] = {
        1.05,  1.0,    0.0,   0.0,   0.0,   0.0,    0.0,   0.0,   0.0,  0.0,
        0.0075, -0.015, 0.0,   0.0,   0.0,   0.0,    0.0,   0.0,   0.0,  0.0,
        -0.9,  0.0,    1.0,   1.0,   0.0,   0.0,    0.0,   0.0,   0.0,  0.0,
        0.0,   0.0,    0.017, -0.04, 0.0,   0.0,    0.0,   0.0,   0.0,  0.0,
        0.0,   0.0,    -0.85, 0.0,   1.0,   1.0,    0.0,   0.0,   0.0,  0.0,
        0.0,   0.0,    0.0,   0.0,   0.04,  -0.103, 0.0,   0.0,   0.0,  0.0,
        0.0,   0.0,    0.0,   0.0,   -0.92, 0.0,    1.0,   1.0,   0.0,  0.0,
        0.0,   0.0,    0.0,   0.0,   0.0,   0.0,    0.105, -0.25, 0.0,  0.0,
        0.0,   0.0,    0.0,   0.0,   0.0,   0.0,    -0.88, 0.0,   1.0,  1.0,
        0.0,   0.0,    0.0,   0.0,   0.0,   0.0,    0.0,   0.0,   0.25, -0.808};
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            (*A)(i, j) = a[i * 10 + j];

    // 2. 10×1 右端向量
    auto b = std::make_unique<Vector>(
        Vector{505.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    // 3. 打印系数矩阵
    A->print();
    // Create SystemEquation with default SORSolver
    SystemEquation systemEq;

    // Set the equation
    systemEq.setEquation(std::move(A), std::move(b));

    // create and set a LUSolver
    auto systemSolver 
                    = std::make_unique<LUSolver>();
    systemEq.setSolver(std::move(systemSolver));

    // Solve the equation
    systemEq.solve();

    // Get and print the solution
    const Vector& solution = systemEq.getSolution();
    std::cout << "Solution:" << std::endl;
    for (const auto& val : solution)
    {
        std::cout << std::fixed << std::setprecision(8) << val << std::endl;
    }
}

int main()
{
    test_system_equation();
    return 0;
}