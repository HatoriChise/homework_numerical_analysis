// SORSolver.cpp

#include <cmath>
#include <iostream>
#include "solver/SORSolver.h"
namespace tinyLinAlg
{

SORSolver::SORSolver(double omega, double tolerance, int maxIterations)
    : omega_(omega), tolerance_(tolerance), maxIterations_(maxIterations)
{
    if(omega_ <= 0.0 || omega_ >= 2.0)
    {
        throw std::invalid_argument("Relaxation factor omega must be in (0, 2)");
    }

    if(tolerance_ <= 0.0)
    {
        throw std::invalid_argument("Tolerance must be positive");
    }

    if(maxIterations_ <= 0)
    {
        throw std::invalid_argument("Maximum number of iterations must be positive");
    }
}

void SORSolver::solve(const Matrix &A, const Vector &b, Vector &x)
{
    // Check Inputs
    if(A.isSquare() == false)
    {
        throw std::invalid_argument("Matrix A must be square.");
    }

    if(A.getRows() != b.size())
    {
        throw std::invalid_argument("Dimension mismatch between A and b.");
    }

    for(size_t i = 0; i < A.getRows(); ++i)
    {
        if(std::abs(A(i, i)) < 1e-15)
        {
            throw std::invalid_argument("Matrix A has zero on diagonal, cannot proceed with SOR.");
        }
    }

    x.resize(A.getRows(), 0.0); // Initialize solution vector x to zero

    for(int iter = 0; iter < maxIterations_; ++iter)
    {
        double maxError = 0.0;

        for(size_t i = 0; i < A.getRows(); ++i)
        {
            double sigma = 0.0;
            for(size_t j = 0; j < A.getCols(); ++j)
            {
                if(j != i)
                {
                    sigma += A(i, j) * x[j];
                }
            }

            double xOld = x[i];
            x[i] = (1 - omega_) * xOld + (omega_ / A(i, i)) * (b[i] - sigma);

            maxError = std::max(maxError, std::abs(x[i] - xOld));
        }

        if(maxError < tolerance_)
        {
            std::cout << "Converged in " << iter + 1 << " iterations." << std::endl;
            return;
        }
    }

    std::cout << "Failed to converge after " << maxIterations_ << " iterations." << std::endl;
    std::cout << "Reached maximum iterations without convergence." << std::endl;
}

} // namespace tinyLinAlg