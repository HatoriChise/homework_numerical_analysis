#include "SystemEquation.h"
#include "Matrix.h"
#include <iostream>

void test_matrix()
{
    tinyLinAlg::Matrix A(3, 3);
    A(0, 0) = 3; A(0, 1) = 2; A(0, 2) = -4;
    A(1, 0) = 2; A(1, 1) = 3; A(1, 2) = 3;
    A(2, 0) = 5; A(2, 1) = -3; A(2, 2) = 1;

    tinyLinAlg::Matrix B(3, 1);
    B(0, 0) = 3;
    B(1, 0) = 15;
    B(2, 0) = 14;

    A.print();
    B.print();

   
}

int main()
{
    test_matrix();
    return 0;
}