#ifndef SYSTEM_EQUATION_H
#define SYSTEM_EQUATION_H

#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "Matrix.h"
#include "solver/LinearSolverBase.h"

namespace tinyLinAlg
{

typedef std::vector<double> Vector;

class SystemEquation
{
private:
    std::unique_ptr<Matrix> matrix_;   // coefficient matrix
    std::unique_ptr<Vector> rhs_;      // right-hand side vector
    std::unique_ptr<Vector> solution_; // solution vector

    std::unique_ptr<LinearSolverBase> solver_;

public:
    explicit SystemEquation(std::unique_ptr<LinearSolverBase> solver = nullptr);

    void setEquation(std::unique_ptr<Matrix> matrix,
                     std::unique_ptr<Vector> rhs);

    void setSolver(std::unique_ptr<LinearSolverBase> solver);
    void solve();

    [[nodiscard]]  Vector& getSolution() const;
    [[nodiscard]]  bool isEquationSet() const
    {
        // if both matrix and rhs are set
        return matrix_ != nullptr && rhs_ != nullptr;
    }
};
} // namespace tinyLinAlg

#endif // SYSTEM_EQUATION_H