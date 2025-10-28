#include "Matrix.h"

namespace tinyLinAlg
{

Matrix::Matrix() : data_(), rows_(0), cols_(0)
{
}

Matrix::Matrix(size_t rows, size_t cols, double initialValue)
    : data_(rows, std::vector<double>(cols, initialValue)), rows_(rows), cols_(cols)
{
}

const double& Matrix::operator()(size_t row, size_t col) const
{
    if (row >= rows_ || col >= cols_)
    {
        throw std::out_of_range("Matrix indices are out of range.");
    }
    return data_.at(row).at(col);
}

double& Matrix::operator()(size_t row, size_t col)
{
    return const_cast<double&>(static_cast<const Matrix&>(*this)(row, col));
}

void Matrix::resize(size_t newRows, size_t newCols, double initialValue)
{
    data_.resize(newRows);
    for (auto& row : data_)
    {
        row.resize(newCols, initialValue);
    }
    rows_ = newRows;
    cols_ = newCols;
}

std::vector<double> Matrix::flatten() const
{
    std::vector<double> flattened;
    flattened.resize(rows_ * cols_);

    for (const auto& row : data_)
    {
        flattened.insert(flattened.end(), row.begin(), row.end());
    }
    return flattened;
} 

void Matrix::print() const
{
    for (const auto& row : data_)
    {
        std::cout << "| ";
        for (const auto& elem : row)
        {
            std::cout << elem << " ";
        }
        std::cout << "| ";
        std::cout << std::endl;
    }
}

}
