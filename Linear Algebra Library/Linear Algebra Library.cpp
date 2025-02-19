// Linear Algebra Library.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "framework.h"
#include "LinearAlgebra.h"

void printMemoryUsage() {
    PROCESS_MEMORY_COUNTERS memInfo;
    GetProcessMemoryInfo(GetCurrentProcess(), &memInfo, sizeof(memInfo));
    std::cout << " Max memory used: " << memInfo.PeakWorkingSetSize / (1024 * 1024) << " MB" << std::endl;
}


double operator*(const Vector& v1, const Vector& v2) {
    double result = 0.0;
	for (int i = 0; i < v1.size(); ++i) {
		result += v1[i] * v2[i];
	}
    return result; 
}


Vector zeros(int n) {
    return Vector(n, 0.0); 
}

// Creates a matrix filled with zeros
Matrix zeros(int rows, int cols) {
    return Matrix(rows, Vector(cols, 0.0));
}

// This function creates a vector full of ones.
Vector ones(int n) {
    return Vector(n, 1.0);
}

// Creates a matrix filled with ones
Matrix ones(int rows, int cols) {
    return Matrix(rows, Vector(cols, 1.0));
}

// This function creates a tensor full of ones.
Tensor ones(int rows, int cols, int n) {
    return Tensor(rows, Matrix(cols, Vector(n, 1.0)));
}

// This function creates a tesseract full of ones.
Tesseract ones(int rows, int cols, int p, int q) {
    return Tesseract(rows, Tensor(cols, Matrix(p, Vector(q, 1.0))));
}

// Creates a column vector of ones of given size
Matrix column_vector(int size) {
    return Matrix(size, Vector(1, 1.0));
}

// Creates an identity matrix of given size
Matrix identity(int size) {
    Matrix I(size, Vector(size, 0.0));
    for (int i = 0; i < size; ++i) {
        I[i][i] = 1.0;
    }
    return I;
}

// This function creates a vector full of random variables.
Vector random(int n) { 
    Vector result(n, 0.0);

    for (int j = 0; j < n; ++j) {
        result[j] = rand(); 
    }
    return result;
}

// This function creates a matrix full of random variables.
Matrix random(int rows, int cols) {
    Matrix result(rows, Vector(cols, 0.0));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i][j] = rand();
        }
    }

    return result;
}

// This function creates a tensor full of random variables.
Tensor random(int rows, int cols, int n) { 
    Tensor result(rows, Matrix(cols, Vector(n, 0.0))); 

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < n; ++k) {
                result[i][j][k] = rand(); 
            }
        }
    }

    return result; 
}

// This function creates a tensor full of random variables.
Tesseract random(int rows, int cols, int p, int q) {
    Tesseract result(rows, Tensor(cols, Matrix(p, Vector(q, 0.0)))); 

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < p; ++k) {
                for (int l = 0; l < q; ++l) {
                    result[i][j][k][l] = rand(); 
                }
            }
        }
    }

    return result;
}

void displayVector(const Vector& A) {
        for (int j = 0; j < A.size(); ++j) {
            cout << A[j] << " "; 
        }
}

// Print normal Matrices;
void displayMatrix(const Matrix& A) {
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[0].size(); ++j) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

bool NanorInf(const Tensor& U) { 
    for (int i = 0; i < U.size(); ++i) { 
        for (int j = 0; j < U[0].size(); ++j) { 
            for (int k = 0; k < U[0][0].size(); ++k) {  
                if (isnan(U[i][j][k]) || isinf(U[i][j][k])) {                  
                    return true; 
                }                
            }
        }
    }    
    return false; 
}

////// Matrix - Matrix LU division

// Performs LU decomposition: A = LU
void LUDecomposition(const Matrix& A, Matrix& L, Matrix& U) {
    int n = A.size();
    L.assign(n, Vector(n, 0.0));
    U.assign(n, Vector(n, 0.0));

    for (int i = 0; i < n; i++) {
        L[i][i] = 1.0; // Diagonal of L is always 1

        // Compute U
        for (int j = i; j < n; j++) {
            U[i][j] = A[i][j];
            for (int k = 0; k < i; k++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }

        // Compute L
        for (int j = i + 1; j < n; j++) {
            if (U[i][i] == 0.0) {
                throw std::runtime_error("Singular matrix detected! LU decomposition failed.");
            }
            L[j][i] = A[j][i];
            for (int k = 0; k < i; k++) {
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
    }
}

// Solves Ly = B using forward substitution
Matrix forwardSubstitution(const Matrix& L, const Matrix& B) {
    int n = L.size();
    int m = B[0].size(); // Number of right-hand sides
    Matrix Y(n, Vector(m, 0.0));

    for (int i = 0; i < n; i++) {
        for (int col = 0; col < m; col++) {
            Y[i][col] = B[i][col];
            for (int j = 0; j < i; j++) {
                Y[i][col] -= L[i][j] * Y[j][col];
            }
        }
    }
    return Y;
}

// Solves Ux = Y using backward substitution
Matrix backwardSubstitution(const Matrix& U, const Matrix& Y) {
    int n = U.size();
    int m = Y[0].size();
    Matrix X(n, Vector(m, 0.0));

    for (int i = n - 1; i >= 0; i--) {
        for (int col = 0; col < m; col++) {
            if (U[i][i] == 0.0) {
                throw std::runtime_error("Singular matrix detected in back-substitution.");
            }
            X[i][col] = Y[i][col];
            for (int j = i + 1; j < n; j++) {
                X[i][col] -= U[i][j] * X[j][col];
            }
            X[i][col] /= U[i][i];
        }
    }
    return X;
}

// Solves AX = B using LU decomposition (overloaded operator)
Matrix operator/(const Matrix& B, const Matrix& A) {
    int n = A.size();
    if (B.size() != n) {
        throw std::invalid_argument("Matrix dimensions do not match for solving AX = B.");
    }

    Matrix L, U;
    LUDecomposition(A, L, U);
    Matrix Y = forwardSubstitution(L, B);
    return backwardSubstitution(U, Y);
}


///// Matrix - Vector LU division

// Solves Ly = B using forward substitution
Vector forwardSubstitution(const Matrix& L, const Vector& B) {
    int n = L.size();
    Vector Y(n, 0.0);

    for (int i = 0; i < n; i++) {
        Y[i] = B[i];
        for (int j = 0; j < i; j++) {
            Y[i] -= L[i][j] * Y[j];
        }
    }
    return Y;
}

// Solves Ux = Y using backward substitution
Vector backwardSubstitution(const Matrix& U, const Vector& Y) {
    int n = U.size();
    Vector X(n, 0.0);

    for (int i = n - 1; i >= 0; i--) {
        if (U[i][i] == 0.0) {
            throw std::runtime_error("Singular matrix detected in back-substitution.");
        }
        X[i] = Y[i];
        for (int j = i + 1; j < n; j++) {
            X[i] -= U[i][j] * X[j];
        }
        X[i] /= U[i][i];
    }
    return X;
}

// Solves AX = B using LU decomposition (overloaded operator)
Vector operator/(const Vector& B, const Matrix& A) { 
    int n = A.size();
    if (B.size() != n) {
        throw std::invalid_argument("Matrix and vector dimensions do not match for solving Ax = B.");
    }
    Matrix L, U;
    LUDecomposition(A, L, U);
    Vector Y = forwardSubstitution(L, B);
    return backwardSubstitution(U, Y);
}

Tensor operator+(const Tensor& A, const Tensor& B) {

    Tensor result(A.size(), Matrix(A[0].size(), Vector(A[0][0].size()))); 

    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[0].size(); ++j) {
            for (int k = 0; k < A[0][0].size(); ++k) {
                result[i][j][k] = A[i][j][k] + B[i][j][k];
            }
        }
    }

    return result; 
}

// Defines operator for matrix multiplication
Matrix operator*(const Matrix& A, const Matrix& B) {
    int Arows = A.size();
    int Acols = A[0].size();
    int Brows = B.size();
    int Bcols = B[0].size();

    if (Acols != Brows) {
        throw std::invalid_argument("Matrix dimensions are not compatible for multiplication.");
    }

    Matrix result(Arows, vector<double>(Bcols, 0.0));
    for (int i = 0; i < Arows; ++i) {
        for (int j = 0; j < Bcols; ++j) {
            for (int k = 0; k < Acols; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

// Defines operator for matrix addition
Matrix operator+(const Matrix& A, const Matrix& B) {
    int Arows = A.size();
    int Acols = A[0].size();
    int Brows = B.size();
    int Bcols = B[0].size();

    if (Acols != Bcols || Arows != Brows) {
        throw std::invalid_argument("Matrix dimensions are not compatible for addition.");
    }

    Matrix result(Arows, Vector(Bcols, 0.0));

    for (int i = 0; i < Arows; ++i) {
        for (int j = 0; j < Bcols; ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }

    return result;
}

// Defines operator for matrix subtraction
Matrix operator-(const Matrix& A, const Matrix& B) {
    int Arows = A.size();
    int Acols = A[0].size();
    int Brows = B.size();
    int Bcols = B[0].size();

    if (Acols != Bcols || Arows != Brows) {
        throw std::invalid_argument("Matrix dimensions are not compatible for subtraction.");
    }

    // Initialize result matrix with proper size
    Matrix result(Arows, Vector(Bcols, 0.0));

    for (int i = 0; i < Arows; ++i) {
        for (int j = 0; j < Bcols; ++j) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }

    return result;
}

// Defines operator for multiplying a matrix by a scalar
Matrix operator*(const double& s, const Matrix& A) {
    int rows = A.size();
    int cols = A[0].size();
    Matrix result(rows, vector<double>(cols, 0));

    // Multiply each element of the matrix by the scalar
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i][j] = s * A[i][j]; // Perform scalar multiplication
        }
    }

    return result;
}

// Defines reverse multiplication of matrix by scalar
Matrix operator*(const Matrix& B, const double& s) {
    return s * B;
}

Vector operator*(const Vector& A, const double& s) {
    int n = A.size(); 
    Vector result(n, 0.0); 
  
    for (int i = 0; i < n; ++i) { 
            result[i] = s * A[i]; // Perform scalar multiplication
    }

    return result;
}

Vector operator*(const double& s, const Vector& A) {
    return A * s; 
}

Matrix operator^(const Matrix& A, const int s) {
    Matrix result(A.size(), Vector(A[0].size()));
    for (int k = 0; k < s; ++k) {
        for (int i = 0; i < A.size(); ++i) {
            for (int j = 0; j < A[0].size(); ++j) { 
                result[i][j] = result[i][j] * result[i][j];  
            }
        }    
    }    
    return result; 
}

Vector operator*(const Matrix& A, const Vector& B) {
    int Arows = A.size();
    int Acols = A[0].size();
    int Brows = B.size();
    int Bcols = 1;

    if (Acols != Brows) {
        throw std::invalid_argument("Matrix dimensions are not compatible for multiplication.");
    }

    Vector result(B.size());

    for (int i = 0; i < Arows; ++i) {
        for (int k = 0; k < Acols; ++k) {
            result[i] += A[i][k] * B[k];
        }
    }

    return result;
}

Vector operator+(const Vector& v1, const Vector& v2) {
    if (v1.size() != v2.size()) {
        throw invalid_argument("Vectors must be the same size");
    }

    Vector result(v1.size());
    for (size_t i = 0; i < v1.size(); i++) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

Vector operator-(const Vector& v1, const Vector& v2) {
    if (v1.size() != v2.size()) { 
        throw invalid_argument("Vectors must be the same size"); 
    } 

    Vector result(v1.size()); 
    for (size_t i = 0; i < v1.size(); i++) { 
        result[i] = v1[i] - v2[i]; 
    }
    return result; 
}

Vector operator/(const Vector& v1, const double& s) {
    Vector result(v1.size()); 
    for (size_t i = 0; i < v1.size(); i++) { 
        result[i] = v1[i]/s ;
    }

    return result;
}


///////////////////////////////////////////////////////////////////////////////////////////
// Matrix1D Class Functions
// Constructor to initialize the matrix with dimensions (k, i)
Matrix1D::Matrix1D(int n, int Nx) : n(n), Nx(Nx) {
    U.resize(n, Vector(Nx, 0.0));
}

// Print Matrix1D
void Matrix1D::print() const {

    for (int i = 0; i < Nx; ++i) {
        cout << "[ ";
        for (int k = 0; k < n; ++k) {
            cout << U[k][i] << " ";
        }
        cout << "] ";
    }
    cout << endl;

}

// This function gets the number of variables
int Matrix1D::nvars() const {
    return U.size();
}

// This function gets the number of columns
int Matrix1D::cols() const {
    return U[0].size();
}

// VariablesProxy constructor
Matrix1D::VariablesProxy::VariablesProxy(Matrix& U, int columnIndex) : U(U), columnIndex(columnIndex) {}

// VariablesProxy assignment operator
Matrix1D::VariablesProxy& Matrix1D::VariablesProxy::operator=(const Matrix& colToSet) {
    for (int i = 0; i < U.size(); ++i) {
        U[i][columnIndex] = colToSet[i][0];
    }
    return *this;
}

// VariablesProxy conversion operator 
Matrix1D::VariablesProxy::operator Matrix() const {
    Matrix column(U.size(), Vector(1, 0.0));
    for (int i = 0; i < U.size(); ++i) {
        column[i][0] = U[i][columnIndex];
    }
    return column;
}

// VariablesProxy subtraction operator
Matrix Matrix1D::VariablesProxy::operator-(const Matrix& other) const {

    Matrix column(U.size(), Vector(1, 0.0));
    for (int i = 0; i < U.size(); ++i) {
        column[i][0] = U[i][columnIndex] - other[i][0]; // Subtract each entry
    }
    return column;
}

// VariabelsProxy .vars caller
Matrix1D::VariablesProxy Matrix1D::vars(int col) {
    return VariablesProxy(U, col);
}

// This function allows the += operator to be used between Matrix1D variables
Matrix1D& Matrix1D::operator+=(const Matrix1D& other) {
    for (int i = 0; i < other.nvars(); ++i)
        for (int j = 0; j < Nx; ++j)
            U[i][j] += other[i][j];

    return *this;
}


// Matrix2D Class Functions // 
// Constructor to initialize the matrix with dimensions (k, i, j)
Matrix2D::Matrix2D(int n, int Nx, int Ny) : n(n), Nx(Nx), Ny(Ny) {
    U.resize(n, Matrix(Nx, Vector(Ny, 0.0)));
}

// Access individual matrix entry at (k, i, j)
double& Matrix2D::entry(int k, int i, int j) {
    return U[k][i][j];
}

// Double checks that the entire matrix does not have infinite or NaN values
bool Matrix2D::NanorInf() const {
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < n; ++k) {
                if (isnan(U[k][i][j]) || isinf(U[k][i][j])) {
                    return true;  // Found NaN or Inf, return true
                    cout << "Inf or NaN found!";
                }
            }
        }
    }
    return false;  // No NaN or Inf found, return false
}

// This function print all elements of Matrix2D
void Matrix2D::print() const {

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            cout << "[ ";
            for (int k = 0; k < n; ++k) {
                cout << U[k][i][j] << " ";
            }
            cout << "] ";
        }
        cout << endl;
    }
}

// This function finds the number of variables
int Matrix2D::nvars() const {
    return U.size();
}

// This function finds the number of rows
int Matrix2D::rows() const {
    return U[0].size();
}

// This function finds the number of columns
int Matrix2D::columns() const {
    return U[0][0].size();
}

// This function allows the += to be used between two variables of type Matrix2D (mainly for U + dU)
Matrix2D& Matrix2D::operator+=(const Matrix2D& other) {
    for (int i = 0; i < other.rows(); ++i) {
        for (int j = 0; j < other.columns(); ++j) {
            for (int k = 0; k < other.nvars(); ++k) {
                U[k][i][j] += other[k][i][j];
            }
        }
    }

    return *this;
}

// This function is the VariablesProxy constructor
Matrix2D::VariablesProxy::VariablesProxy(Tensor& U, int i, int j) : U(U), i(i), j(j) {}

// This function converts VariablesProxy class to a Matrix class (get all variables at (i, j))
Matrix2D::VariablesProxy::operator Matrix() const {
    Matrix result(U.size(), Vector(1, 0.0));
    for (int k = 0; k < U.size(); ++k) {
        result[k][0] = U[k][i][j];  // Extract all variable values     
    }
    return result;
}

// This function allows setting assigning matrices to the Matrix2D ".vars" function.
Matrix2D::VariablesProxy& Matrix2D::VariablesProxy::operator=(const Matrix& values) {
    for (int k = 0; k < U.size(); ++k) {
        U[k][i][j] = values[k][0];  // Assign new values        
    }
    return *this;
}

// This function create the ".vars" caller
Matrix2D::VariablesProxy Matrix2D::vars(int i, int j) {
    return VariablesProxy(U, i, j);
}

// This function is the LineProxy constructor
Matrix2D::LineProxy::LineProxy(Tensor& U, int i) : U(U), i(i) {};

// This function converts the ".line" function into a Matrix if necessary
Matrix2D::LineProxy::operator Matrix() const {
    Matrix result(U.size(), Vector(U[0][0].size(), 0.0));
    for (int j = 0; j < U[0][0].size(); ++j) {
        for (int k = 0; k < U.size(); ++k) {
            result[k][j] = U[k][i][j];
        }
    }
    return result;
}

// This function allows setting a ".line" equal to a matrix
Matrix2D::LineProxy& Matrix2D::LineProxy::operator=(const Matrix& values) {
    for (int j = 0; j < U[0][0].size(); ++j) {
        for (int k = 0; k < U.size(); ++k) {
            U[k][i][j] = values[k][j];
        }
    }
    return *this;
}

// This function allows setting a ".line" function equal to a Matrix1D class.
Matrix2D::LineProxy& Matrix2D::LineProxy::operator=(const Matrix1D& values) {
    for (int k = 0; k < U.size(); ++k) {
        for (int j = 0; j < U[0][0].size(); ++j) {
            U[k][i][j] = values[k][j];
        }
    }

    return*this;
}

// LineProxy .line caller
Matrix2D::LineProxy Matrix2D::line(int i) {
    return LineProxy(U, i);
}


// This function is the StripProxy constructor
Matrix2D::StripProxy::StripProxy(Tensor& U, int j) : U(U), j(j) {};

// This function converts a ".strip" function into a matrix
Matrix2D::StripProxy::operator Matrix() const {
    Matrix result(U.size(), Vector(U[0].size(), 0.0));
    for (int i = 0; i < U[0].size(); ++i) {
        for (int k = 0; k < U.size(); ++k) {
            result[k][i] = U[k][i][j];
        }
    }
    return result;
}

// This function allows setting a ".strip" function equal to a matrix
Matrix2D::StripProxy& Matrix2D::StripProxy::operator=(const Matrix& values) {
    for (int i = 0; i < U[0].size(); ++i) {
        for (int k = 0; k < U.size(); ++k) {
            U[k][i][j] = values[k][i];
        }
    }
    return *this;
}

// StripProxy .line caller
Matrix2D::StripProxy Matrix2D::strip(int j) {
    return StripProxy(U, j);
}