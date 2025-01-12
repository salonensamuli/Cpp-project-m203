// Paris Dauphine-PSL M203 introduction to c++ 
// course home assignment by Samuli Salonen (22301060)
// main source: https://rtraba.com/wp-content/uploads/2015/05/solving_pdes_in_c.pdf

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
using namespace std;

class Matrix { // (chapter 2)
private:
    vector<vector<double>> data_;
    int rows_, cols_;

public:
    // constructors
    Matrix(int rows, int cols); // to initialize a matrix
    Matrix(const Matrix& other); // to copy matrices

    // accessors
    int rows() const;
    int cols() const;

    // operator() to access matrix elements
    double& operator()(int i, int j); // to modify matrix elements
    const double& operator()(int i, int j) const; // to declare matrix objects

    // utility functions
    bool isSquare() const; // to check if the matrix is square
    double determinant() const; // to calculate the determinant (chapter 2)
    bool isInvertible() const; // to check if the matrix is invertible
    Matrix inverse() const; // to calcule the inverse matrix (chapter 2)
};

#endif