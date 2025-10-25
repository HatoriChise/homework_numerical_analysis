#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include <vector>

class CubicSpline
{
public:
    // Constructor
    CubicSpline(const std::vector<double>& x, const std::vector<double>& y);

    // Destructor
    ~CubicSpline();

    // Main solving function
    void solve();

    // Evaluate spline at given x value
    double evaluate(double x) const;

    // Get spline coefficients for debugging/analysis
    std::vector<std::vector<double>> getCoefficients() const;

    // Set boundary conditions
    void setBoundaryConditions(double left_slope,
                               double right_slope); // Clamped
    void setNaturalBoundary(); // Natural (second derivative = 0 at endpoints)

private:
    std::vector<double> xValues_; // x coordinates of data points
    std::vector<double> yValues_; // y coordinates of data points
    std::vector<double> hValues_; // intervals (h_i = x_{i+1} - x_i)

    // Spline coefficients for each interval: S_i(x) = a_i + b_i*(x-x_i) +
    // c_i*(x-x_i)^2 + d_i*(x-x_i)^3
    std::vector<double> aCoeffs_; // a coefficients
    std::vector<double> bCoeffs_; // b coefficients
    std::vector<double> cCoeffs_; // c coefficients
    std::vector<double> dCoeffs_; // d coefficients

    int numberOfPoints_;          // number of data points
    bool naturalBoundary_;        // flag for natural boundary condition
    double leftSlope_, rightSlope_; // boundary slopes for clamped conditions
    bool boundarySet_;            // flag indicating if boundary conditions were set

    // Helper methods
    void setupSystem();      // Set up the tridiagonal system
    void solveTridiagonal(); // Solve the tridiagonal system for c coefficients
    void computeCoefficients(); // Compute b and d coefficients from c
                                // coefficients
};

#endif // CUBIC_SPLINE_H