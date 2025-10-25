#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cubic_spline.h"

int main()
{
    std::cout << "Cubic Spline Interpolation Implementation" << std::endl;

    // Example usage of cubic spline
    // Define sample data points (x, y)
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> y = {0.0, 1.0, 0.0, -1.0, 0.0};

    // Create cubic spline object
    CubicSpline spline(x, y);

    // Solve for spline coefficients
    spline.solve();

    // Example: Evaluate spline at specific points
    std::cout << "\nSpline evaluation:" << std::endl;
    for (double t = 0.0; t <= 4.0; t += 0.5)
    {
        double value = spline.evaluate(t);
        std::cout << "f(" << t << ") = " << std::fixed << std::setprecision(6)
                  << value << std::endl;
    }

    return 0;
}