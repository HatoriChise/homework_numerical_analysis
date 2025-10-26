#include <algorithm>
#include <cmath>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include "cubic_spline.h"

CubicSpline::CubicSpline(const std::vector<double>& x,
                         const std::vector<double>& y)
    : xValues_(std::move(x)), yValues_(std::move(y)),
      numberOfPoints_(xValues_.size())
{
    // Validate inputs (size, finiteness, strictly increasing x)
    checkInputs();

    // Initialize coefficient vectors
    aCoeffs_.resize(numberOfPoints_ - 1);
    bCoeffs_.resize(numberOfPoints_ - 1);
    cCoeffs_.resize(numberOfPoints_ - 1);
    dCoeffs_.resize(numberOfPoints_ - 1);

    hValues_.resize(numberOfPoints_ - 1);

    // Calculate h values (intervals) — increasing order already validated
    for (std::size_t i = 0; i + 1 < static_cast<std::size_t>(numberOfPoints_);
         ++i)
    {
        hValues_[i] = xValues_[i + 1] - xValues_[i];
    }

    // Initialize boundary conditions
    naturalBoundary_ = true;
    boundarySet_ = false;
}

void CubicSpline::solve()
{
    // Set up the system of equations based on boundary conditions
    setupSystem();

    // Solve the tridiagonal system to find c coefficients
    // TDMA algorithm
    solveTridiagonal();

    // Compute b and d coefficients
    computeCoefficients();
}

void CubicSpline::setupSystem()
{
    const auto n = numberOfPoints_;

    switch (boundaryType_)
    {
    case BoundaryType::Natural:
        // Natural boundary conditions (second derivatives at endpoints are
        // zero)
        setupSystemNatural();
        break;

    case BoundaryType::Clamped:
        // Clamped boundary conditions (first derivatives at endpoints are
        // specified)
        setupSystemClamped();
        break;

    case BoundaryType::NotAKnot:
        // Not-a-knot boundary conditions (third derivative is continuous at the
        // second point)
        setupSystemNotAKnot();
        break;

    default:
        throw std::invalid_argument("Invalid boundary condition");
        break;
    }
}

void CubicSpline::solveTridiagonal()
{
    // Solve the tridiagonal system using the Thomas algorithm
    // After solving, c_coeffs will contain the c coefficients
}

void CubicSpline::computeCoefficients()
{
    // Compute a, b, and d coefficients from the c coefficients
    // a[i] = y[i]
    // b[i] = (y[i+1] - y[i])/h[i] - h[i]*(2*c[i] + c[i+1])/3
    // d[i] = (c[i+1] - c[i])/(3*h[i])
}

double CubicSpline::evaluate(double x) const
{
    // Find the appropriate interval [x[i], x[i+1]] that contains x
    // Then evaluate the cubic polynomial S_i(x) = a[i] + b[i]*(x-x[i]) +
    // c[i]*(x-x[i])^2 + d[i]*(x-x[i])^3
    return 0.0; // Placeholder
}

std::vector<std::vector<double>> CubicSpline::getCoefficients() const
{
    std::vector<std::vector<double>> result(4);
    result[0] = aCoeffs_;
    result[1] = bCoeffs_;
    result[2] = cCoeffs_;
    result[3] = dCoeffs_;
    return result;
}

void CubicSpline::setBoundaryConditions(BoundaryType type, double firstValue,
                                        double secondValue)
{
    boundaryType_ = type;
    boundarySet_ = true;

    switch (type)
    {
    case BoundaryType::Natural:
        // 不需要额外参数
        naturalBoundary_ = true;
        break;

    case BoundaryType::Clamped:
        naturalBoundary_ = false;
        leftSlope_ = firstValue;
        rightSlope_ = secondValue;
        break;

    case BoundaryType::NotAKnot:
        // 不需要额外参数
        naturalBoundary_ = false;
        break;
    }
}


void CubicSpline::checkInputs() const
{
    const auto n = xValues_.size();

    if (n != yValues_.size())
    {
        throw std::invalid_argument(
            "Invalid input: x and y must have the same size");
    }
    if (n < 2)
    {
        throw std::invalid_argument(
            "Invalid input: need at least 2 data points");
    }

    // Check finiteness (no NaN/Inf)
    for (auto x : xValues_)
    {
        if (!std::isfinite(x))
        {
            throw std::invalid_argument("Invalid input: x must be finite");
        }
    }
    for (auto y : yValues_)
    {
        if (!std::isfinite(y))
        {
            throw std::invalid_argument("Invalid input: y must be finite");
        }
    }

    // Check strictly increasing x values
    for (std::size_t i = 1; i < n; ++i)
    {
        if (xValues_[i + 1] <= xValues_[i])
        {
            std::ostringstream oss;
            oss << "x[" << i << "] = " << xValues_[i] << " is not less than x["
                << i + 1 << "] = " << xValues_[i + 1];
            throw std::invalid_argument(oss.str());
        }
    }
}

void CubicSpline::setupSystemNatural()
{

}

void CubicSpline::setupSystemClamped()
{

}

void CubicSpline::setupSystemNotAKnot()
{
    
}