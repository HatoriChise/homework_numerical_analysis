#ifndef LINEAR_SOLVER_BASE_H
#define LINEAR_SOLVER_BASE_H

#include <vector>
#include "Matrix.h"

namespace tinyLinAlg
{

using Vector = std::vector<double>;


/**
 * @brief  Abstract base class for linear solvers.
 *  Defines the interface for solving linear systems of equations of the form Ax = b.
 */
class LinearSolverBase
{
public:

    virtual ~LinearSolverBase() = default;

    /**
     * @brief Solve the linear system Ax = b.
     *
     * @param A The matrix A in Ax = b.
     * @param b The vector b in Ax = b.
     * @param x The solution vector x.
     */
    virtual void solve(const Matrix &A, const Vector &b, Vector &x) = 0;
};

} // namespace tinyLinAlg

#endif // LINEAR_SOLVER_BASE_H