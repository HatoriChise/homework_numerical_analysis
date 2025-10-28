#include "Matrix.h"

namespace tinyLinAlg
{

template <typename T> Matrix<T>::Matrix() : rows_(0), cols_(0)
{
}

template <typename T>
Matrix<T>::Matrix(size_t rows, size_t cols, T initialValue)
    : rows_(rows), cols_(cols), data_(rows, std::vector<T>(cols, initialValue))
{
}

template <typename T>
const T& Matrix<T>::operator()(size_t row, size_t col) const
{
    if (row >= rows_ || col >= cols_)
    {
        throw std::out_of_range("Matrix indices are out of range.");
    }
    return data_.at(row).at(col);
}

template <typename T>
 T& Matrix<T>::operator()(size_t row, size_t col)
{
    return const_cast<T&>(static_cast<const Matrix<T>&>(*this)(row, col));
}

template <typename T>
void Matrix<T>::resize(size_t newRows, size_t newCols, T initialValue)
{
    data_.resize(newRows);
    for (auto& row : data_)
    {
        row.resize(newCols, initialValue);
    }
    rows_ = newRows;
    cols_ = newCols;
}

template <typename T>
std::vector<T> Matrix<T>::flatten() const
{
    std::vector<T> flattened;
    flattened.resize(rows_ * cols_);

    for (const auto& row : data_)
    {
        flattened.insert(flattened.end(), row.begin(), row.end());
    }
    return flattened;
} 

template <typename T>
void Matrix<T>::print() const
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

}// namespace tinyLinAlg

