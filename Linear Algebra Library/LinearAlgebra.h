#pragma once

#include <iostream>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <stdexcept>  // For exception handling
#include <Windows.h>
#include <Psapi.h>
#include <chrono>
#include <unordered_map>
#include <string>
#include <array>

using namespace std;
constexpr double pi = 3.141592653; 

void printMemoryUsage(); 

typedef array<double, 4> Vector;   
typedef array<Vector, 4> Matrix; 
  
Vector zerosV();
Matrix zerosM();
Vector onesV();
Matrix onesM(); 
Matrix identity();
Vector randomV();
Matrix randomM();

void displayVector(const Vector& A);
void displayMatrix(const Matrix& A);

double operator*(const Vector& A, const Vector& B);  
Vector operator+(const Vector& v1, const Vector& v2);
Vector operator-(const Vector& v1, const Vector& v2); 
Vector operator*(const Vector& A, const double& s);
Vector operator*(const double& s, const Vector& A); 
Vector operator/(const Vector& v1, const double& s);


Matrix outerProduct(const Vector& a, const Vector& b);

Vector operator*(const Matrix& A, const Vector& B);


Matrix operator+(const Matrix& A, const Matrix& B);
Matrix operator-(const Matrix& A, const Matrix& B);
Matrix operator*(const Matrix& A, const double& s);
Matrix operator*(const double& s, const Matrix& A); 
Matrix operator*(const Matrix& A, const Matrix& B);



void LUDecomposition(const Matrix& A, Matrix& L, Matrix& U);
Vector forwardSubstitution(const Matrix& L, const Vector& B);
Vector backwardSubstitution(const Matrix& U, const Vector& Y); 
Vector operator/(const Vector& B, const Matrix& A); 


Matrix forwardSubstitution(const Matrix& L, const Matrix& b);
Matrix backwardSubstitution(const Matrix& U, const Matrix& y);
Matrix operator/(const Matrix& b, const Matrix& A);

