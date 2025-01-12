// Paris Dauphine-PSL M203 introduction to c++ 
// course home assignment by Samuli Salonen (22301060)
// main source: https://rtraba.com/wp-content/uploads/2015/05/solving_pdes_in_c.pdf

#include "matrix.h"
#include <stdexcept> // for errors
#include <cmath> // for abs
#include <algorithm> // for swap
using namespace std;

// constructors, see .h for further explanations
Matrix::Matrix(int rows, int cols)
    : rows_(rows), cols_(cols), data_(rows, vector<double>(cols, 0.0)) {}
Matrix::Matrix(const Matrix& other)
    : rows_(other.rows_), cols_(other.cols_), data_(other.data_) {}

// accessors
int Matrix::rows() const {
    return rows_;
}
int Matrix::cols() const {
    return cols_;
}

// () to access matrix elements, see .h for further explanations
double& Matrix::operator()(int i, int j) {
    return data_[i][j];
}
const double& Matrix::operator()(int i, int j) const {
    return data_[i][j];
}

// check if square
bool Matrix::isSquare() const {
    return (rows_ == cols_);
}

// calculate the determinant
double Matrix::determinant() const {
    // check if square
    if (!isSquare()) {
        throw runtime_error("Matrix not square! (determinant)");
    }

    Matrix A(*this); // copy the original matrix
    double det = 1.0; // initialize determinant
    int n = rows_; // get number of rows and columns

    // Gaussian elimination (chapter 2)
    for (int i = 0; i < n; ++i) {
        double pivot = A(i, i); // pivot element (diagonal in current column)
        if (abs(pivot) < 1e-14) { // if the pivot is near zero (unstable), we need to swap rows to find a non-zero pivot
            bool pivot_found = false;

            // search for a non-zero pivot in rows below the current row
            for (int k = i + 1; k < n; ++k) {
                if (abs(A(k, i)) > 1e-14) {
                    swap(A.data_[i], A.data_[k]); // swap the current row with the row having a valid pivot
                    det = -det; // change the determinant's sign due to the row swap
                    pivot_found = true;
                    break;
                }
            }
            pivot = A(i, i); // update the pivot after the swap

            // if no valid pivot was found, the matrix is singular (determinant is zero)
            if (!pivot_found || abs(pivot) < 1e-14) {
                return 0.0; // singular
            }
        }

        det *= pivot; // multiply the determinant by the pivot value

        // normalize the current row by dividing by the pivot
        for (int r = i + 1; r < n; ++r) {
            double factor = A(r, i) / pivot;
            // subtract the scaled row from rows below to eliminate the current column
            for (int c = i; c < n; ++c) {
                A(r, c) -= factor * A(i, c);
            }
        }
    }
    return det;
}


// check if invertible
bool Matrix::isInvertible() const {
    if(!isSquare()) return false; // false if square
    double det = determinant();
    return (abs(det) > 1e-14); // false if singular
}

// calculate the inverse
Matrix Matrix::inverse() const {
    // check if square
    if (!isSquare()) {
        throw runtime_error("Matrix not square! (inverse)");
    }

    int n = rows_; // number of rows and columns
    Matrix A(*this); // copy the original matrix
    Matrix I(n, n);  // identity matrix
    for (int i = 0; i < n; ++i) {
        I(i, i) = 1.0; // diagonal elemetns to 1
    }

    // Gaussian elimination (again chapter 2): transforme A into identity matrix and I into the inverse of the original
    for (int i = 0; i < n; ++i) {
        double pivot = A(i, i); // pivot element (diagonal in current column)
        if (abs(pivot) < 1e-14) { // if the pivot is near zero (unstable), swap rows to find a non-zero pivot
            bool pivot_found = false;

            // search for a non-zero pivot in rows below the current row
            for (int k = i + 1; k < n; ++k) {
                if (abs(A(k, i)) > 1e-14) {
                    // swap the current row with the row having a valid pivot
                    swap(A.data_[i], A.data_[k]);
                    swap(I.data_[i], I.data_[k]);
                    pivot_found = true;
                    break; // stop if valid pivot
                }
            }
            pivot = A(i, i); // update the pivot after the swap

            // if no valid pivot was found, the matrix is singular (determinant is zero)
            if (!pivot_found || abs(pivot) < 1e-14) {
                throw runtime_error("Singular matrix, no inverse!");
            }
        }

        // normalize the current row by dividing by the pivot
        for (int c = 0; c < n; ++c) {
            A(i, c) /= pivot;
            I(i, c) /= pivot;
        }
        // eliminate all other entries in the current column (make them zero)
        for (int r = 0; r < n; ++r) {
            if (r != i) {
                double factor = A(r, i);
                // subtract the scaled row from rows below to eliminate the current column
                for (int c = 0; c < n; ++c) {
                    A(r, c) -= factor * A(i, c);
                    I(r, c) -= factor * I(i, c);
                }
            }
        }
    }
    return I; // now A is an identity matrix and I should be the inverse :)
}