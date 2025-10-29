// LUSolver.cpp
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "solver/LUSolver.h"

namespace tinyLinAlg
{

LUSolver::LUSolver() : isDecomposed_(false)
{
}

LUSolver::LUSolver(const Matrix &A) : isDecomposed_(false)
{
    bool success = decompose(A);
    if(!success)
    {
        throw std::runtime_error("LU decomposition failed during construction.");
    }
}

void LUSolver::printLU() const
{
    std::cout << "L Matrix:" << std::endl;
    LMatrix_.print();
    std::cout << "U Matrix:" << std::endl;
    UMatrix_.print();

}
bool LUSolver::decompose(const Matrix &A)
{
    std::size_t n = A.getRows();
    // Perform LU decomposition
    LMatrix_ = Matrix(n, n, 0.0);
    UMatrix_ = Matrix(n, n, 0.0);

    for (auto i =0; i < n; ++i)
    {
        for(auto j = i; j < n; ++j)
        {
            auto sum = 0.0;
            // Upper Triangular
            for (auto k = 0; k < i; ++k)
            {
                sum += LMatrix_(i, k) * UMatrix_(k, j);
            }
            UMatrix_(i, j) = A(i, j) - sum;

        }
        LMatrix_(i, i) = 1.0; // Diagonal as 1

        for (auto j = i + 1; j < n; ++j)
        {
            auto sum = 0.0;
            // Lower Triangular
            for (auto k = 0; k < i; ++k)
            {
                sum += LMatrix_(j, k) * UMatrix_(k, i);
            }
            
            LMatrix_(j, i) = (A(j, i) - sum) / UMatrix_(i, i);
        }
    }

    isDecomposed_ = true; // Set the flag to indicate that decomposition is successful
    printLU(); // Print the L and U matrices
    return true;
}

void LUSolver::forwardSubstitution(const Vector &b, Vector &y)
{
    y = Vector(b.size(), 0.0); // Initialize y with zeros
    for (std::size_t i = 0; i < b.size(); ++i)
    {
        double sum = 0.0;
        for (std::size_t j = 0; j < i; ++j)
        {
            sum += LMatrix_(i, j) * y[j];
        }
        y[i] = b[i] - sum;
    }
}

void LUSolver::backwardSubstitution(const Vector &y, Vector &x)
{
    x = Vector(y.size(), 0.0); // Initialize x with zeros
    for (std::size_t i = y.size(); i-- > 0;)
    {
        double sum = 0.0;
        for (std::size_t j = i + 1; j < y.size(); ++j)
        {
            sum += UMatrix_(i, j) * x[j];
        }
        if (UMatrix_(i, i) == 0.0)
        {
            throw std::runtime_error("Zero pivot encountered during backward substitution.");
        }
        x[i] = (y[i] - sum) / UMatrix_(i, i);
    }
}

void LUSolver::solve(const Matrix &A, const Vector &b, Vector &x)
{
    if(!decompose(A))
    {
        throw std::runtime_error("LU decomposition failed.");
    }
    Vector y;
    forwardSubstitution(b, y);
    backwardSubstitution(y, x);
}



} // namespace tinyLinAlg
