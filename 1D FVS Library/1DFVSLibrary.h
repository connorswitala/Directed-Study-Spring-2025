#pragma once


#include <iostream>
#include <vector>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream> 

#include "LinearAlgebra.h"


constexpr double gamma = 1.4;

typedef pair<Matrix, Matrix> Duo;

double calculate_dt(Matrix1D& U, const double& dx, const int& Nx, const double CFL);

Matrix primtoCons(const Matrix& V);
Matrix constoPrim(const Matrix& U);

Matrix A_Minus(const Matrix& U, int nx);
Matrix A_Plus(const Matrix& U, int nx);

Duo Boundary1D(const Matrix U, int nx, const string type);
Matrix1D solveLine(Matrix1D U, const double dx, const int Nx, const double dt, Matrix A, Matrix B, Matrix C, Matrix F, Matrix alpha, Matrix2D v, Matrix2D g);

void saveVectorsToCSV(const string& filename,
	const Vector& density,
	const Vector& u_vel,
	const Vector& pressure,
	const Vector& x);

