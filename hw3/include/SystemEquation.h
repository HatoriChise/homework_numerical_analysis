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
    std::unique_ptr<Matrix> matrix_; // coefficient matrix
    std::unique_ptr<Vector> rhs_; // right-hand side vector 
    std::unique_ptr<Vector> solution_; // solution vector
public:

};
} // namespace tinyLinAlg

#endif // SYSTEM_EQUATION_H