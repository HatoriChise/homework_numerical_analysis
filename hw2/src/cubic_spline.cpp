#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include "cubic_spline.h"

CubicSpline::CubicSpline(const std::vector<double>& x,
                         const std::vector<double>& y)
    : xValues_(x), yValues_(y), hValues_(), aCoeffs_(), bCoeffs_(), cCoeffs_(),
      dCoeffs_(), lowerDiag_(), mainDiag_(), upperDiag_(), rhs_(),
      numberOfPoints_(static_cast<int>(x.size())), boundarySet_(false),
      naturalBoundary_(true), boundaryType_(BoundaryType::Natural)
{
    // Validate inputs (size, finiteness, strictly increasing x)
    checkInputs();

    // Initialize coefficient vectors
    // eg. numberOfPoints_ = 6, aCoeffs_ = [a0, a1, a2, a3, a4, a5] (6 elements)
    aCoeffs_.resize(numberOfPoints_);
    bCoeffs_.resize(numberOfPoints_ - 1);
    cCoeffs_.resize(numberOfPoints_);
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
    solution_ = solveTridiagonal();

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

std::vector<double> CubicSpline::solveTridiagonal()
{
    // Solve the tridiagonal system using the Thomas algorithm
    // After solving, c_coeffs will contain the c coefficients
    const auto systemSize = mainDiag_.size();


    std::vector<double> cPrime(systemSize);
    std::vector<double> dPrime(systemSize);
    std::vector<double> x(systemSize); // solution vector

    // Forward sweep
    // step 1: modify the first row (i=0)
    auto denmo = mainDiag_[0];
    if (std::abs(denmo) < 1e-10)
    {
        throw std::runtime_error("Division by zero in Thomas algorithm");
    }
    cPrime[0] = upperDiag_[0] / denmo;
    dPrime[0] = rhs_[0] / denmo;

    // step 2: modify the remaining rows
    for (auto i = 1; i < systemSize - 1; ++i)
    {
        denmo = mainDiag_[i] - lowerDiag_[i - 1] * cPrime[i - 1];
        if (std::abs(denmo) < 1e-10)
        {
            throw std::runtime_error("Division by zero in Thomas algorithm");
        }
        cPrime[i] = upperDiag_[i] / denmo;
        dPrime[i] = (rhs_[i] - lowerDiag_[i - 1] * dPrime[i - 1]) / denmo;
    }

    // step 3: last row (i = systemSize - 1)
    denmo = mainDiag_[systemSize - 1] -
            lowerDiag_[systemSize - 2] * cPrime[systemSize - 2];
    if (std::abs(denmo) < 1e-10)
    {
        throw std::runtime_error("Division by zero in Thomas algorithm");
    }
    dPrime[systemSize - 1] =
        (rhs_[systemSize - 1] -
         lowerDiag_[systemSize - 2] * dPrime[systemSize - 2]) /
        denmo;

    // Backward sweep
    x[systemSize - 1] = dPrime[systemSize - 1];
    for (auto i = systemSize - 2; i >= 1; --i)
    {
        x[i] = dPrime[i] - cPrime[i] * x[i + 1];
    }
    x[0] = dPrime[0] - cPrime[0] * x[1];

    return x;
}

void CubicSpline::computeCoefficients()
{
    // Compute the coefficients a, b, c, d for each spline segment

    aCoeffs_ = yValues_; // a_i = y_i
    if (boundaryType_ == BoundaryType::Natural)
    {
        cCoeffs_[0] = 0.0;
        cCoeffs_[numberOfPoints_ - 1] = 0.0;

        for (std::size_t i = 0; i < solution_.size(); ++i)
        {
            cCoeffs_[i + 1] = solution_[i] / 2;
        }
    }
    else if (boundaryType_ == BoundaryType::Clamped)
    {
        cCoeffs_ = solution_;
    }

    for (std::size_t i = 0; i < numberOfPoints_ - 1; ++i)
    {
        std::cout << "c_" << i << " = " << cCoeffs_[i] << std::endl;
        std::cout << "y_" << i << " = " << yValues_[i] << std::endl;
        std::cout << "h_" << i << " = " << hValues_[i] << std::endl;

        bCoeffs_[i] = (yValues_[i + 1] - yValues_[i]) / hValues_[i]
                    - (hValues_[i]/3) * (cCoeffs_[i + 1] + 2*cCoeffs_[i]);
        std::cout << "b_" << i << " = " << bCoeffs_[i] << std::endl;
        dCoeffs_[i] = (cCoeffs_[i + 1] - cCoeffs_[i]) / (3.0 * hValues_[i]);
    }
}


double CubicSpline::evaluate(double x) const
{
    if (x < xValues_[0] || x > xValues_[numberOfPoints_ - 1]) {
        throw std::out_of_range("Input x is outside the range of the interpolation points.");
    }


    auto it = std::upper_bound(xValues_.begin(), xValues_.end(), x);
    std::size_t i = std::distance(xValues_.begin(), it) - 1;
    
    if (i >= numberOfPoints_ - 1) {
        i = numberOfPoints_ - 2; // 确保在最后一个区间
    }
    
    double dx = x - xValues_[i];
    double dx2 = dx * dx;
    double dx3 = dx2 * dx;

    // 4. 计算并返回 S_i(x)
    return aCoeffs_[i] + bCoeffs_[i] * dx + cCoeffs_[i] * dx2 + dCoeffs_[i] * dx3;
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
        const auto fBackward =
            (yValues_[i] - yValues_[i - 1]) / hValues_[i - 1];
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
    // todo: Set up the system of equations for not-a-knot boundary conditions
    const auto n = numberOfPoints_;

    hValues_.resize(n - 1);
    for (std::size_t i = 0; i < n - 1; ++i)
    {
        hValues_[i] = xValues_[i + 1] - xValues_[i];
    }
}

void CubicSpline::printDebugInfo() const
{
    // Lambda to print a vector with its name and size
    auto print_vec = [](const char* name, const std::vector<double>& v) {
        std::cout << name << " (size=" << v.size() << "):";
        if (v.empty())
        {
            std::cout << " []" << std::endl;
            return;
        }
        std::cout << " [";
        for (std::size_t i = 0; i < v.size(); ++i)
        {
            std::cout << std::fixed << std::setprecision(9) << v[i];
            if (i + 1 < v.size())
                std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    };

    // Lambda to convert BoundaryType to string
    auto boundary_to_string = [this]() -> const char* {
        switch (boundaryType_)
        {
        case BoundaryType::Natural:
            return "Natural";
        case BoundaryType::Clamped:
            return "Clamped";
        case BoundaryType::NotAKnot:
            return "NotAKnot";
        default:
            return "Unknown";
        }
    };

    std::cout << "===== CubicSpline Debug Info =====" << std::endl;
    std::cout << "numberOfPoints_: " << numberOfPoints_ << std::endl;
    std::cout << "boundarySet_: " << (boundarySet_ ? "true" : "false")
              << ", naturalBoundary_: " << (naturalBoundary_ ? "true" : "false")
              << std::endl;
    std::cout << "boundaryType_: " << boundary_to_string() << std::endl;
    std::cout << "leftSlope_: " << leftSlope_
              << ", rightSlope_: " << rightSlope_ << std::endl;

    print_vec("xValues_", xValues_);
    print_vec("yValues_", yValues_);
    print_vec("hValues_", hValues_);

    // Tridiagonal system
    print_vec("mainDiag_", mainDiag_);
    print_vec("lowerDiag_", lowerDiag_);
    print_vec("upperDiag_", upperDiag_);
    print_vec("rhs_", rhs_);
    // solution vector (c coefficients)
    print_vec("solution_ (c coefficients)", solution_);
    // Coefficients

    print_vec("aCoeffs_", aCoeffs_);
    print_vec("bCoeffs_", bCoeffs_);
    print_vec("cCoeffs_", cCoeffs_);
    print_vec("dCoeffs_", dCoeffs_);

    std::cout << "==================================" << std::endl;
}