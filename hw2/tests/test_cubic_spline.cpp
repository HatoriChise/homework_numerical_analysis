#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "../include/cubic_spline.h"

namespace
{

/**
 * @brief Expects an invalid_argument to be thrown by the given function.
 *
 * @param name The name of the test.
 * @param fn The function to test.
 * @return true If the function throws an invalid_argument.
 * @return false Otherwise.
 */
bool expectThrowsInvalidArgument(const std::string& name,
                                 const std::function<void()>& fn)
{
    try
    {
        fn();
        std::cout << "  FAIL: " << name
                  << " — expected invalid_argument, but no exception thrown"
                  << std::endl;
        return false;
    }
    catch (const std::invalid_argument& e)
    {
        std::cout << "  PASS: " << name
                  << " — caught invalid_argument: " << e.what() << std::endl;
        return true;
    }
    catch (const std::exception& e)
    {
        std::cout << "  FAIL: " << name
                  << " — caught unexpected exception: " << e.what()
                  << std::endl;
        return false;
    }
}

bool expectNoThrow(const std::string& name, const std::function<void()>& fn)
{
    try
    {
        fn();
        std::cout << "  PASS: " << name << std::endl;
        return true;
    }
    catch (const std::exception& e)
    {
        std::cout << "  FAIL: " << name
                  << " — unexpected exception: " << e.what() << std::endl;
        return false;
    }
}

} // namespace

int main()
{
    std::cout << "Running CubicSpline input validation tests" << std::endl;
    std::cout << "==========================================" << std::endl;

    int passed = 0;
    int total = 0;

    // 1) Valid minimal input
    ++total;
    passed += expectNoThrow("valid minimal input (2 points)", [] {
        std::vector<double> x{0.0, 1.0};
        std::vector<double> y{0.0, 1.0};
        CubicSpline s(x, y);
    });

    // 2) Size mismatch
    ++total;
    passed += expectThrowsInvalidArgument("size mismatch", [] {
        std::vector<double> x{0.0, 1.0, 2.0};
        std::vector<double> y{0.0, 1.0};
        CubicSpline s(x, y);
    });

    // 3) Too few points
    ++total;
    passed += expectThrowsInvalidArgument("too few points", [] {
        std::vector<double> x{0.0};
        std::vector<double> y{0.0};
        CubicSpline s(x, y);
    });

    // 4) Non-increasing x (equal)
    ++total;
    passed += expectThrowsInvalidArgument("non-increasing (equal)", [] {
        std::vector<double> x{0.0, 1.0, 1.0};
        std::vector<double> y{0.0, 1.0, 2.0};
        CubicSpline s(x, y);
    });

    // 5) Non-increasing x (decreasing)
    ++total;
    passed += expectThrowsInvalidArgument("non-increasing (decreasing)", [] {
        std::vector<double> x{0.0, 2.0, 1.0};
        std::vector<double> y{0.0, 4.0, 1.0};
        CubicSpline s(x, y);
    });

    // 6) NaN in x
    ++total;
    passed += expectThrowsInvalidArgument("NaN in x", [] {
        std::vector<double> x{0.0, std::numeric_limits<double>::quiet_NaN(),
                              1.0};
        std::vector<double> y{0.0, 1.0, 2.0};
        CubicSpline s(x, y);
    });

    // 7) Inf in x
    ++total;
    passed += expectThrowsInvalidArgument("Inf in x", [] {
        std::vector<double> x{0.0, std::numeric_limits<double>::infinity(),
                              1.0};
        std::vector<double> y{0.0, 1.0, 2.0};
        CubicSpline s(x, y);
    });

    // 8) NaN in y
    ++total;
    passed += expectThrowsInvalidArgument("NaN in y", [] {
        std::vector<double> x{0.0, 1.0, 2.0};
        std::vector<double> y{0.0, std::numeric_limits<double>::quiet_NaN(),
                              2.0};
        CubicSpline s(x, y);
    });

    // 9) Inf in y
    ++total;
    passed += expectThrowsInvalidArgument("Inf in y", [] {
        std::vector<double> x{0.0, 1.0, 2.0};
        std::vector<double> y{0.0, std::numeric_limits<double>::infinity(),
                              2.0};
        CubicSpline s(x, y);
    });

    // 10) Boundary API usable (Natural default, and Clamped with finite slopes)
    ++total;
    passed += expectNoThrow("boundary API calls (no solve)", [] {
        std::vector<double> x{0.0, 1.0, 2.0};
        std::vector<double> y{0.0, 1.0, 4.0};
        CubicSpline s(x, y);
        s.setBoundaryConditions(CubicSpline::BoundaryType::Natural, 0.0, 0.0);
        s.setBoundaryConditions(CubicSpline::BoundaryType::Clamped, 1.0, 2.0);
        s.setBoundaryConditions(CubicSpline::BoundaryType::NotAKnot, 0.0, 0.0);
    });

    std::cout << "\nTests passed: " << passed << "/" << total << std::endl;

    return (passed == total) ? 0 : 1;
}
