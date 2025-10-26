#include <algorithm>
#include <cmath>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include "cubic_spline.h"

CubicSpline::CubicSpline(const std::vector<double>& x,
                         const std::vector<double>& y)
    : xValues_(x), yValues_(y), hValues_(), aCoeffs_(), bCoeffs_(), cCoeffs_(),
      dCoeffs_(), lowerDiag_(), mainDiag_(), upperDiag_(), rhs_(),
      numberOfPoints_(static_cast<int>(x.size())), boundarySet_(false),
      naturalBoundary_(true)
{
    // Default boundary: Natural
    boundaryType_ = BoundaryType::Natural;
    leftSlope_ = rightSlope_ = 0.0;
    leftSecondDerivative_ = rightSecondDerivative_ = 0.0;

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

    // Check strictly increasing x values using adjacent comparison
    for (std::size_t i = 0; i + 1 < n; ++i)
    {
        if (!(xValues_[i + 1] > xValues_[i]))
        {
            std::ostringstream oss;
            oss << "Invalid input: x must be strictly increasing. Found x[" << i
                << "] = " << xValues_[i] << ", x[" << (i + 1)
                << "] = " << xValues_[i + 1];
            throw std::invalid_argument(oss.str());
        }
    }
}

void CubicSpline::setupSystemClamped()
{
    // Set up the system of equations for clamped boundary conditions
    const auto n = numberOfPoints_;

    hValues_.resize(n - 1);
    for (std::size_t i = 0; i < n - 1; ++i)
    {
        hValues_[i] = xValues_[i + 1] - xValues_[i];
    }

    const std::size_t systemSize = n;
    mainDiag_.resize(systemSize);
    lowerDiag_.resize(systemSize - 1);
    upperDiag_.resize(systemSize - 1);
    rhs_.resize(systemSize);

    // first boundary equation (i=0)
    mainDiag_[0] = 2.0 * hValues_[0];
    upperDiag_[0] = hValues_[0];
    rhs_[0] = 6.0 * ((yValues_[1] - yValues_[0]) / hValues_[0] - leftSlope_);
    
    // last boundary equation (i=n-1)
    mainDiag_[n - 1] = 2.0 * hValues_[n - 2];
    lowerDiag_[n - 2] = hValues_[n - 2];
    rhs_[n - 1] = 6.0 * (rightSlope_ -
                         (yValues_[n - 1] - yValues_[n - 2]) / hValues_[n - 2]);

    for (std::size_t i = 1; i < systemSize - 1; ++i)
    {
        mainDiag_[i] = 2.0 * (hValues_[i - 1] + hValues_[i]);
        lowerDiag_[i - 1] = hValues_[i - 1];
        upperDiag_[i] = hValues_[i];

        const auto fForward = (yValues_[i + 1] - yValues_[i]) / hValues_[i];
        const auto fBackward = (yValues_[i] - yValues_[i - 1]) / hValues_[i - 1];
        rhs_[i] = 6.0 * (fForward - fBackward);
    }
}

void CubicSpline::setupSystemNatural()
{
    // 在自然边界条件中，M0 = 0 和 Mn = 0，只需求解中间的 n-2 个方程
    // Set up the system of equations for natural boundary conditions
    const auto n = numberOfPoints_; // we have n data points

    hValues_.resize(n - 1); 
    for (std::size_t i = 0; i < n - 1; ++i)
    {
        hValues_[i] = xValues_[i + 1] - xValues_[i];
    }

    // allocate space for the tridiagonal matrix and right-hand side
    const std::size_t systemSize = n - 2;
    mainDiag_.resize(systemSize);
    lowerDiag_.resize(systemSize - 1);
    upperDiag_.resize(systemSize - 1);
    rhs_.resize(systemSize);

    for (std::size_t i = 0; i < systemSize; ++i)
    {
        // Fill in the tridiagonal matrix and rhs based on natural conditions
        // 因为我们只解中间的 n-2 个方程，索引需要偏移 1
        // i 对应方程中索引 i 的方程，对应数据点索引 dataIndex = i + 1
        const std::size_t dataIndex = i + 1;
        // main diagonal = 2 * (h[i] + h[i+1])
        mainDiag_[i] = 2.0 * (hValues_[dataIndex - 1] + hValues_[dataIndex]);
        // rhs = 6 * (y[i+1] - y[i]) / h[i] - 6 * (y[i] - y[i-1]) / h[i-1]
        const auto fForward = (yValues_[dataIndex + 1] - yValues_[dataIndex]) /
                              hValues_[dataIndex];
        const auto fBackward = (yValues_[dataIndex] - yValues_[dataIndex - 1]) /
                               hValues_[dataIndex - 1];
        rhs_[i] = 6.0 * (fForward - fBackward);

        // upper and lower diagonals
        if (i < systemSize - 1)
        {
            upperDiag_[i] = hValues_[dataIndex];
        }
        if (i > 0)
        {
            lowerDiag_[i - 1] = hValues_[dataIndex - 1];
        }
    }
}


void CubicSpline::setupSystemNotAKnot()
{
    const auto n = numberOfPoints_;

    hValues_.resize(n - 1);
    for (std::size_t i = 0; i < n - 1; ++i)
    {
        hValues_[i] = xValues_[i + 1] - xValues_[i];
    }
}