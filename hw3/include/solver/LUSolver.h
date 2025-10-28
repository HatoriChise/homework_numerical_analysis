#ifndef LU_SOLVER_H
#define LU_SOLVER_H

#include "LinearSolverBase.h"

namespace tinyLinAlg
{
class LUSolver : public LinearSolverBase
{
public:
    void solve(const Matrix &A, const Vector &b, Vector &x) override;
};
} // namespace tinyLinAlg

#endif