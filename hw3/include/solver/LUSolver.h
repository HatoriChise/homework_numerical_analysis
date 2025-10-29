#ifndef LU_SOLVER_H
#define LU_SOLVER_H

#include "LinearSolverBase.h"

namespace tinyLinAlg
{
class LUSolver : public LinearSolverBase
{

private:
    Matrix LMatrix_;
    Matrix UMatrix_;

    bool isDecomposed_ = false; // flag to check if the matrix has been decomposed

    bool decompose(const Matrix &A);
    void forwardSubstitution(const Vector &b, Vector &y);
    void backwardSubstitution(const Vector &y, Vector &x);
public:

    LUSolver() ;

    explicit LUSolver(const Matrix &A);

    /**
     * @brief  Solve the system Ax = b.
     *
     * @param A Coefficient matrix.
     * @param b Right-hand side vector.
     * @param x Solution vector.
     */
    void solve(const Matrix &A, const Vector &b, Vector &x) override;

    void printLU() const;
};
} // namespace tinyLinAlg

#endif