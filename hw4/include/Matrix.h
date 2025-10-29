#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>
#include <iostream>

namespace tinyLinAlg
{

class Matrix
{ 
private:

    std::vector<std::vector<double>> data_;
    size_t rows_;
    size_t cols_;

public:
    // default constructor
    Matrix() ;

    // constructor with initial value, rows and columns
    Matrix(size_t rows, size_t cols, double initialValue = 0.0);

    virtual ~Matrix() = default;

    // get number of rows and columns
    std::size_t getRows() const
    {
        return rows_;
    }

    std::size_t getCols() const
    {
        return cols_;
    }

    bool isSquare() const
    {
        return rows_ == cols_;
    }


    // get element at (row, col)
    const double& operator()(size_t row, size_t col) const;
    double& operator()(size_t row, size_t col);
    void resize(size_t newRows, size_t newCols, double initialValue = 0.0);
    std::vector<double> flatten() const;
    void print() const;
};



} // namespace tinyLinAlg
#endif // MATRIX_H
