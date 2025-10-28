#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>

namespace tinyLinAlg
{

template <typename T>
class Matrix
{ 
private:
    using Matrix = std::vector<std::vector<T>>;
    Matrix data_;
    size_t rows_;
    size_t cols_;

public:
    // default constructor
    Matrix() = default;

    // constructor with initial value, rows and columns
    Matrix(size_t rows, size_t cols, T initialValue = T())
    {
        data_.resize(rows, std::vector<T>(cols, initialValue));
    }

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
    const T& operator()(size_t row, size_t col) const;
    T& operator()(size_t row, size_t col);

    void resize(size_t newRows, size_t newCols, T initialValue = T());
    std::vector<T> flatten() const;

    void print() const;
};

} // namespace tinyLinAlg
#endif // MATRIX_H
