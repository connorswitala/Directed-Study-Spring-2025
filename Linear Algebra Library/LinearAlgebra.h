#pragma once

#include <iostream>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <stdexcept>  // For exception handling

using namespace std;

typedef vector<vector<vector<double>>> Tensor;
typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

Matrix operator*(const Matrix& A, const Matrix& B);
Matrix operator+(const Matrix& A, const Matrix& B);
Matrix operator-(const Matrix& A, const Matrix& B);
Matrix operator*(const double& s, const Matrix& B);
Matrix operator*(const Matrix& B, const double& s);
Matrix operator^(const Matrix& A, const int s);   

void displayMatrix(const Matrix& A);

Matrix zeros(int rows, int cols);
Matrix ones(int rows, int cols);
Matrix identity(int size);
Matrix column_vector(int size);
Matrix random(int rows, int cols);


void LUDecomposition(const Matrix& A, Matrix& L, Matrix& U);
Matrix forwardSubstitution(const Matrix& L, const Matrix& b);
Matrix backwardSubstitution(const Matrix& U, const Matrix& y);
Matrix operator/(const Matrix& b, const Matrix& A);

class Matrix1D {
private:
    Matrix U;
    int n, Nx;

public:
    Matrix1D(int n, int Nx);
    void print() const;
    int cols() const;
    int nvars() const;

    Vector& operator[](int i) { return U[i]; }
    const Vector& operator[](int i) const { return U[i]; }

    class VariablesProxy {
    private:

        Matrix& U;
        int columnIndex;

    public:
        VariablesProxy(Matrix& U, int columnIndex);
        VariablesProxy& operator=(const Matrix& colToSet);
        operator Matrix() const;
        Matrix operator-(const Matrix& other) const;

    };

    VariablesProxy vars(int col);
    Matrix1D& operator+=(const Matrix1D& other);

};

class Matrix2D {
private:
    Tensor U;
    int n, Nx, Ny;

public:

    Matrix2D(int n, int Nx, int Ny);
    double& entry(int k, int i, int j);
    bool NanorInf() const;
    void print() const;
    int nvars() const;
    int rows() const;
    int columns() const;

    Matrix& operator[](int i) { return U[i]; }
    const Matrix& operator[](int i) const { return U[i]; }

    class VariablesProxy {
    private:
        Tensor& U;
        int i, j;

    public:

        VariablesProxy(Tensor& U, int i, int j);

        // Convert ColumnProxy to a Vector (get all variables at (i, j))
        operator Matrix() const;

        // Assign a Vector to all variables at (i, j)
        VariablesProxy& operator=(const Matrix& values);
    };


    class LineProxy {
    private:
        Tensor& U;
        int i;

    public:

        LineProxy(Tensor& U, int i);
        operator Matrix() const;
        LineProxy& operator=(const Matrix& values);
        LineProxy& operator=(const Matrix1D& values);
    };

    class StripProxy {
    private:
        Tensor& U;
        int j;

    public:

        StripProxy(Tensor& U, int j);
        StripProxy(const Tensor& U, int j);
        operator Matrix() const;
        StripProxy& operator=(const Matrix& values);

    };

    // Return ColumnProxy for (i, j) that holds all variables
    VariablesProxy vars(int i, int j);
    LineProxy line(int i);
    StripProxy strip(int j);
    Matrix2D& operator+=(const Matrix2D& other);


};
