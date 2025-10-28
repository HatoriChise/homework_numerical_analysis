#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cubic_spline.h"

void test()
{
    std::vector<double> x {0.0, 20.0, 40.0, 60.0, 80.0, 100.0}; 
    std::vector<double> y {3.5, 4.2, 3.8, 4.5, 4.0, 3.6};
    CubicSpline spline(x, y);
    auto boundaryType = CubicSpline::BoundaryType::Natural;
    spline.setBoundaryConditions(boundaryType, 0.0, 0.0);
    spline.solve();
    spline.printDebugInfo();
    // 进行插值测试
    std::cout << "Evaluating spline at various points:" << std::endl;
    for (double xi = 0.0; xi <= 100.0; xi += 10.0)
    {
        std::cout << "x = " << xi << " y = " << spline.evaluate(xi) << std::endl;
    }
}

int main()
{
    test();
    return 0;
}