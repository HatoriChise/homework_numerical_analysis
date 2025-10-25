#include <cmath>
#include <stdexcept>
#include "cubic_spline.h"

CubicSpline::CubicSpline(const std::vector<double>& x,
                         const std::vector<double>& y)
{
    if (x.size() != y.size() || x.size() < 2)
    {
        throw std::invalid_argument("Invalid input: x and y must have the same "
                                    "size and at least 2 points");
    }

    // Copy input data
    x_values = x;
    y_values = y;
    n = x.size();

    // Initialize coefficient vectors
    a_coeffs.resize(n - 1);
    b_coeffs.resize(n - 1);
    c_coeffs.resize(n);
    d_coeffs.resize(n - 1);

    h_values.resize(n - 1);

    // Calculate h values (intervals)
    for (int i = 0; i < n - 1; ++i)
    {
        h_values[i] = x_values[i + 1] - x_values[i];
        if (h_values[i] <= 0)
        {
            throw std::invalid_argument(
                "x values must be in strictly increasing order");
        }
    }

    // Initialize boundary conditions
    natural_boundary = true;
    boundary_set = false;
}

CubicSpline::~CubicSpline()
{
    // Cleanup if needed
}

void CubicSpline::solve()
{
    // Set up the system of equations based on boundary conditions
    setupSystem();

    // Solve the tridiagonal system to find c coefficients
    solveTridiagonal();

    // Compute b and d coefficients
    computeCoefficients();
}

void CubicSpline::setupSystem()
{
    // Set up the tridiagonal matrix system Ac = r for the c coefficients
    // This method will populate the matrix coefficients based on boundary
    // conditions

    if (natural_boundary)
    {
        // Natural spline: second derivatives at endpoints are zero
        // This means c[0] = 0 and c[n-1] = 0
    }
    else
    {
        // Clamped spline: derivatives at endpoints are given
        // This means we'll have different equations for the first and last rows
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
    result[0] = a_coeffs;
    result[1] = b_coeffs;
    result[2] = c_coeffs;
    result[3] = d_coeffs;
    return result;
}

void CubicSpline::setBoundaryConditions(double left_slope, double right_slope)
{
    natural_boundary = false;
    left_slope = left_slope;
    right_slope = right_slope;
    boundary_set = true;
}

void CubicSpline::setNaturalBoundary()
{
    natural_boundary = true;
    boundary_set = true;
}