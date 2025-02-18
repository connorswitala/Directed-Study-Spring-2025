#pragma once

#include "LinearAlgebra.h"
#include "GridGenerator.h"
#include <cstdlib>
#include <fstream> 
#include <omp.h>

using namespace std;

constexpr double gamma = 1.4;
typedef pair<Vector, Matrix> Duo;

enum class BoundaryCondition {

	Wall, 
	Inlet, 
	Outlet, 
	Symmetry
}; 

struct BoundaryConditions {

	BoundaryCondition left;
	BoundaryCondition right;
	BoundaryCondition bottom;
	BoundaryCondition top;

	BoundaryConditions(BoundaryCondition left, BoundaryCondition right, BoundaryCondition bottom, BoundaryCondition top) : left(left), right(right), bottom(bottom), top(top) {}

};

void writeVTK(const string& filename,
	Tensor& U,
	const Grid& grid,
	const int Nx,
	const int Ny);

void writeCSV(const string& filename,
	const Tensor& U,
	const Grid& grid,
	const int j,
	const int Nx);

Vector primtoCons(const Vector& V);
Vector constoPrim(const Vector& U);


Duo Boundary2D(BoundaryCondition type, const Vector& U, const Vector& U_inlet, const Point& normals);


void solveOneTimestep(Tensor& U, Tensor& dU_new, Vector U_inlet, Tensor& dU_old, const Grid& grid,
	BoundaryConditions BoundaryTypes, const int& Nx, const int& Ny, double& dt, const double& CFL,
	Tensor& i_Fluxes, Tensor& j_Fluxes, Tesseract& i_plus_Jacobians, Tesseract& i_minus_Jacobians, Tesseract& j_plus_Jacobians, Tesseract& j_minus_Jacobians);

double calculate_dt(Tensor& U, const Grid& grid, const int& Nx, const int& Ny, const double& CFL);

void Calculate_Jacobians(const Tensor& U, Vector& U_inlet, const BoundaryConditions& BoundaryTypes, const Grid& grid, int Nx, int Ny,
	Tensor& i_Fluxes, Tensor& j_Fluxes, Tesseract& i_plus_Jacobians, Tesseract& i_minus_Jacobians, Tesseract& j_plus_Jacobians, Tesseract& j_minus_Jacobians); 

void solveLeftLine(Tensor& U, Tensor& dU_new, Vector U_inlet, Tensor& dU_old, const Grid& grid,
	BoundaryConditions BoundaryType, const int i, const int Nx, const int Ny, const double dt,
	Tensor& i_Fluxes, Tensor& j_Fluxes, Tesseract& i_plus_Jacobians, Tesseract& i_minus_Jacobians, Tesseract& j_plus_Jacobians, Tesseract& j_minus_Jacobians);

void solveMiddleLine(Tensor& U, Tensor& dU_new, Vector U_inlet, Tensor& dU_old, const Grid& grid,
	BoundaryConditions BoundaryType, const int i, const int Nx, const int Ny, const double dt,
	Tensor& i_Fluxes, Tensor& j_Fluxes, Tesseract& i_plus_Jacobians, Tesseract& i_minus_Jacobians, Tesseract& j_plus_Jacobians, Tesseract& j_minus_Jacobians);

void solveRightLine(Tensor& U, Tensor& dU_new, Vector U_inlet, Tensor& dU_old, const Grid& grid, 
	BoundaryConditions BoundaryType, const int i, const int Nx, const int Ny, const double dt,
	Tensor& i_Fluxes, Tensor& j_Fluxes, Tesseract& i_plus_Jacobians, Tesseract& i_minus_Jacobians, Tesseract& j_plus_Jacobians, Tesseract& j_minus_Jacobians);


double calculateInnerResidual(Tensor& U, Tensor& dU_new, Tensor& dU_old, const Grid& grid, const int Nx, const int Ny, const double dt,
	Tensor& i_Fluxes, Tensor& j_Fluxes, Tesseract& i_plus_Jacobians, Tesseract& i_minus_Jacobians, Tesseract& j_plus_Jacobians, Tesseract& j_minus_Jacobians);

double calculateResidual(const Tensor& U, const Grid& grid, int Nx, int Ny, const Tensor& i_Fluxes, const Tensor& j_Fluxes); 
