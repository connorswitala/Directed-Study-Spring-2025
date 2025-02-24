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


// HERE //
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



class Solver { 

private:

	CellTensor U, dU_new, dU_old;   

	iFaceTensor i_Fluxes; 
	iFaceTesseract i_plus_Jacobians, i_minus_Jacobians;

	jFaceTensor j_Fluxes;
	jFaceTesseract j_plus_Jacobians, j_minus_Jacobians;

	Grid& grid;   
	double dt, inner_residual; 

	BoundaryConditions BoundaryType; 
	const double CFL;
	Vector V_inlet, U_inlet;  


public: 

	double outer_residual;

	Solver(const Vector& V_inlet, Grid& grid, BoundaryConditions BoundaryType, const double CFL);

	Vector primtoCons(const Vector& V);  
	Vector constoPrim(const Vector& U);

	Matrix Boundary2DE(BoundaryCondition type, const Vector& U, const Point& normals);
	Vector Boundary2DU(BoundaryCondition type, const Vector& U, const Point& normals);

	void solve(); 
	void solveOneTimestep();
	void calculate_dt(); 
	void Calculate_Jacobians(); 

	void solveLeftLine();
	void solveMiddleLine(const int i);
	void solveRightLine(); 

	void calculateInnerResidual();
	void calculateResidual(); 

	void write2DCSV(const string& filename);
	
	void time(void (Solver::* func)()) { 
		auto start = TIME;

		(this->*func)();

		auto end = TIME;
		DURATION duration = end - start;
		cout << "Time taken: " << duration.count() << endl;
	}


};






 