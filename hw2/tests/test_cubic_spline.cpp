#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "../include/cubic_spline.h"

bool testBasicFunctionality()
{
    std::cout << "Testing basic functionality..." << std::endl;

    // Define simple test data - a quadratic function y = x^2
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> y = {0.0, 1.0, 4.0, 9.0}; // y = x^2

    try
    {
        CubicSpline spline(x, y);
        spline.solve();

        // Test a few evaluation points
        std::vector<double> test_points = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
        bool all_pass = true;

        for (double t : test_points)
        {
            double interpolated = spline.evaluate(t);
            double expected = t * t; // x^2
            double error = std::abs(interpolated - expected);

            if (error > 1e-6)
            {
                std::cout << "  FAIL: At x=" << t << ", got " << interpolated
                          << ", expected " << expected << ", error: " << error
                          << std::endl;
                all_pass = false;
            }
        }

        if (all_pass)
        {
            std::cout << "  PASS: Basic functionality test" << std::endl;
            return true;
        }
        else
        {
            std::cout << "  FAIL: Basic functionality test" << std::endl;
            return false;
        }
    }
    catch (const std::exception& e)
    {
        std::cout << "  ERROR: Exception occurred: " << e.what() << std::endl;
        return false;
    }
}

bool testNaturalBoundary()
{
    std::cout << "Testing natural boundary conditions..." << std::endl;

    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {1.0, 2.0, 3.0};

    try
    {
        CubicSpline spline(x, y);
        spline.setNaturalBoundary();
        spline.solve();

        std::cout << "  PASS: Natural boundary conditions test" << std::endl;
        return true;
    }
    catch (const std::exception& e)
    {
        std::cout << "  ERROR: Exception occurred: " << e.what() << std::endl;
        return false;
    }
}

bool testClampedBoundary()
{
    std::cout << "Testing clamped boundary conditions..." << std::endl;

    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 0.0};

    try
    {
        CubicSpline spline(x, y);
        spline.setBoundaryConditions(1.0, -1.0); // Set derivatives at endpoints
        spline.solve();

        std::cout << "  PASS: Clamped boundary conditions test" << std::endl;
        return true;
    }
    catch (const std::exception& e)
    {
        std::cout << "  ERROR: Exception occurred: " << e.what() << std::endl;
        return false;
    }
}

int main()
{
    std::cout << "Running Cubic Spline Tests" << std::endl;
    std::cout << "=========================" << std::endl;

    int passed = 0;
    int total = 3;

    if (testBasicFunctionality())
        passed++;
    std::cout << std::endl;

    if (testNaturalBoundary())
        passed++;
    std::cout << std::endl;

    if (testClampedBoundary())
        passed++;
    std::cout << std::endl;

    std::cout << "Tests passed: " << passed << "/" << total << std::endl;

    if (passed == total)
    {
        std::cout << "All tests PASSED!" << std::endl;
        return 0;
    }
    else
    {
        std::cout << "Some tests FAILED!" << std::endl;
        return 1;
    }
}