#pragma once

#include "LinearAlgebra.h"
#include "GridGenerator.h"
#include <cstdlib>
#include <fstream> 
#include <omp.h>
#include <chrono> 

#define TIME chrono::high_resolution_clock::now(); 
#define DURATION chrono::duration<double> duration; 

using namespace std;

constexpr double gamma = 1.4;
constexpr double R = 287.0;
constexpr double cv = R / (gamma - 1);
constexpr double cp = cv + R;
constexpr double Pr = 0.72;

constexpr int Nx = 100;
constexpr int Ny = 50;



using CellTensor = array < array < Vector, Ny>, Nx>; 
using CellTesseract = array < array < Matrix, Ny>, Nx>;  

using iFaceTensor = array < array < Vector, Ny>, Nx + 1>;  
using iFaceTesseract = array < array <Matrix, Ny>, Nx + 1>; 

using jFaceTensor = array < array < Vector, Ny + 1>, Nx>; 
using jFaceTesseract = array < array <Matrix, Ny + 1>, Nx>;   

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

void write1DCSV(const string& filename,
	const CellTensor& U,
	const Grid& grid,
	const int j); 

void write2DCSV(const string& filename,
	CellTensor& U,
	const Grid& grid); 

Vector primtoCons(const Vector& V);
Vector constoPrim(const Vector& U);


Duo Boundary2D(BoundaryCondition type, const Vector& U, const Vector& U_inlet, const Point& normals);


void solveOneTimestep(CellTensor& U, CellTensor& dU_new, Vector U_inlet, CellTensor& dU_old, const Grid& grid,
	BoundaryConditions BoundaryTypes, double& dt, const double& CFL,
	iFaceTensor& i_Fluxes, jFaceTensor& j_Fluxes, iFaceTesseract& i_plus_Jacobians, iFaceTesseract& i_minus_Jacobians, jFaceTesseract& j_plus_Jacobians, jFaceTesseract& j_minus_Jacobians);

double calculate_dt(CellTensor& U, const Grid& grid, const double& CFL);

void Calculate_Jacobians(const CellTensor& U, Vector& U_inlet, const BoundaryConditions& BoundaryTypes, const Grid& grid,
	iFaceTensor& i_Fluxes, jFaceTensor& j_Fluxes, iFaceTesseract& i_plus_Jacobians, iFaceTesseract& i_minus_Jacobians, jFaceTesseract& j_plus_Jacobians, jFaceTesseract& j_minus_Jacobians);

void solveLeftLine(CellTensor& U, CellTensor& dU_new, Vector U_inlet, CellTensor& dU_old, const Grid& grid,
	BoundaryConditions BoundaryType, const int i, const double dt,
	iFaceTensor& i_Fluxes, jFaceTensor& j_Fluxes, iFaceTesseract& i_plus_Jacobians, iFaceTesseract& i_minus_Jacobians, jFaceTesseract& j_plus_Jacobians, jFaceTesseract& j_minus_Jacobians);

void solveMiddleLine(CellTensor& U, CellTensor& dU_new, Vector U_inlet, CellTensor& dU_old, const Grid& grid,
	BoundaryConditions BoundaryType, const int i, const double dt,
	iFaceTensor& i_Fluxes, jFaceTensor& j_Fluxes, iFaceTesseract& i_plus_Jacobians, iFaceTesseract& i_minus_Jacobians, jFaceTesseract& j_plus_Jacobians, jFaceTesseract& j_minus_Jacobians);

void solveRightLine(CellTensor& U, CellTensor& dU_new, Vector U_inlet, CellTensor& dU_old, const Grid& grid, 
	BoundaryConditions BoundaryType, const int i, const double dt,
	iFaceTensor& i_Fluxes, jFaceTensor& j_Fluxes, iFaceTesseract& i_plus_Jacobians, iFaceTesseract& i_minus_Jacobians, jFaceTesseract& j_plus_Jacobians, jFaceTesseract& j_minus_Jacobians);


double calculateInnerResidual(CellTensor& U, CellTensor& dU_new, CellTensor& dU_old, const Grid& grid, const double dt,
	iFaceTensor& i_Fluxes, jFaceTensor& j_Fluxes, iFaceTesseract& i_plus_Jacobians, iFaceTesseract& i_minus_Jacobians, jFaceTesseract& j_plus_Jacobians, jFaceTesseract& j_minus_Jacobians);

double calculateResidual(const CellTensor& U, const Grid& grid, const iFaceTensor& i_Fluxes, const jFaceTensor& j_Fluxes); 
