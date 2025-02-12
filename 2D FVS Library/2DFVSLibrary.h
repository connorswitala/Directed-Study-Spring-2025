#pragma once

#include "LinearAlgebra.h"
#include "GridGenerator.h"
#include <cstdlib>
#include <fstream> 

using namespace std;

enum class BoundaryCondition {
	Wall, 
	Inlet, 
	Outlet, 
	Symmetry
}; 


constexpr double gamma = 1.4;
typedef pair<Matrix, Matrix> Duo;

double calculate_dt(Matrix2D& U, const Grid& grid, int Nx, int Ny, double CFL);

double calculateInnerResidual(Matrix2D& U, Matrix2D& dU_old, Matrix2D& dU_new, const Grid& grid, const int Nx, const int Ny, const double dt, int i);
double calculateResidual(Matrix2D& U, const Grid& grid, int Nx, int Ny);


Matrix1D solveLeftLine(Matrix2D& U, Matrix U_inlet, Matrix2D& dU_old, const Grid& grid,	BoundaryCondition leftType, BoundaryCondition toptype, BoundaryCondition bottomtype, const int i, const int Nx, const int Ny, const double dt); 
Matrix1D solveMiddleLine(Matrix2D& U, Matrix U_inlet, Matrix2D& dU_old, const Grid& grid, BoundaryCondition toptype, BoundaryCondition bottomtype, const int i, const int Nx, const int Ny, const double dt);
Matrix1D solveRightLine(Matrix2D& U, Matrix U_inlet, Matrix2D& dU_old, const Grid& grid, BoundaryCondition righttype, BoundaryCondition toptype, BoundaryCondition bottomtype, const int i, const int Nx, const int Ny, const double dt);

Matrix2D solveOneTimestep(Matrix2D& U, Matrix2D& dU_old, Matrix U_inlet, const Grid& grid, double dt, int Nx, int Ny, BoundaryCondition left, BoundaryCondition right, BoundaryCondition bottom, BoundaryCondition top); 

Matrix primtoCons(const Matrix& V);
Matrix constoPrim(const Matrix& U);

Matrix A_Minus(const Matrix& U, Point normals);
Matrix A_Plus(const Matrix& U, Point normals);

Duo Boundary2D(BoundaryCondition type, const Matrix U, const Matrix U_inlet, const Point normals); 


void writeVTK(const string& filename,
	Matrix2D& U,
	const Grid& grid,
	const int Nx,
	const int Ny); 

void writeCSV(const string& filename,
	const Matrix& density,
	const Matrix& u_velocity,
	const Matrix& v_velocity,
	const Matrix& pressure,
	const Grid& grid,
	const int n,
	const int Nx);