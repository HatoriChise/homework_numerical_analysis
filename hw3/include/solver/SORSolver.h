#ifndef SOR_SOLVER_H
#define SOR_SOLVER_H

#include "LinearSolverBase.h"

namespace tinyLinAlg
{
class SORSolver : public LinearSolverBase
{
private:
    double omega_; // relaxation factor
    double tolerance_; // convergence tolerance
    int maxIterations_; // maximum number of iterations

public:
    
    /**
     * @brief  Construct a new SORSolver object
     * 
     * @param omega Relaxation factor (0 < omega < 2)
     * @param tolerance Convergence tolerance
     * @param maxIterations Maximum number of iterations
     */
    explicit SORSolver(double omega=1.0, double tolerance=1e-6, int maxIterations=1000);

    /**
     * @brief Solve the system Ax = b using the SOR method
     * 
     * @param A Matrix A coefficient of the system
     * @param b Vector b rhs of the system
     * @param x Vector x solution of the system
     */
    void solve(const Matrix &A, const Vector &b, Vector &x) override;
};
} // namespace tinyLinAlg

#endif // SOR_SOLVER_H