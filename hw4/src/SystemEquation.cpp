#include "SystemEquation.h"
#include "solver/SORSolver.h"
#include "solver/LUSolver.h"

namespace tinyLinAlg
{

SystemEquation::SystemEquation(std::unique_ptr<LinearSolverBase> solver)
    : solver_(std::move(solver))
{
    // If no solver is provided, use a default SORSolver
    if (!solver_)
    {
        solver_ = std::make_unique<SORSolver>(1, 1e-10, 1000);
    }
}

void SystemEquation::setEquation(std::unique_ptr<Matrix> matrix,
                                 std::unique_ptr<Vector> rhs)
{
    if(matrix->isSquare() == false)
    {
        throw std::invalid_argument("Coefficient matrix must be square");
    }
    if (matrix->getRows() != rhs->size())
    {
        throw std::invalid_argument("Matrix and RHS must have the same size");
    }

    this->matrix_ = std::move(matrix);
    this->rhs_ = std::move(rhs);

    solution_.reset();
}

void SystemEquation::setSolver(std::unique_ptr<LinearSolverBase> solver)
{
    // ownership transferred to SystemEquation
    solver_ = std::move(solver);
}
void SystemEquation::solve()
{
    if (!isEquationSet())
    {
        throw std::runtime_error(
            "Equation is not set. Call setEquation() first.");
    }
    if (!solver_)
    {
        throw std::logic_error("Solver is not set. Call setSolver() first.");
    }

    if (!solution_)
    {
        solution_ = std::make_unique<Vector>(matrix_->getRows(), 0.0);
    }
    else if (solution_->size() != rhs_->size())
    {
        throw std::invalid_argument("Solution and RHS must have the same size");
    }

    solver_->solve(*matrix_, *rhs_, *solution_);
    
}
Vector& SystemEquation::getSolution() const
{
    if (!solution_)
    {
        throw std::runtime_error("Solution has not been computed yet.");
    }
    return *solution_;
}

} // namespace tinyLinAlg
