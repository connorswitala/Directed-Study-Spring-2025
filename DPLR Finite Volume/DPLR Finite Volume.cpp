#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <fstream> 
#include "GridGenerator.h" 
#include "2DFVSLibrary.h" 
#include "LinearAlgebra.h"
#include <iomanip>


using namespace std;


#define TIME chrono::high_resolution_clock::now(); 
#define DURATION chrono::duration<double> duration;   

int main() {

	auto start = TIME;

	int counter = 0;
	double CFL = 2.0, dt;

	double p = 10000.0, T = 300.0, R = 287.0, M = 2.5, a = sqrt(gamma * R * T), u = M * a, v = 0, rho = p / (R * T); 

	Vector V_inlet = { rho, u, v, p };
	Vector U_inlet = primtoCons(V_inlet);


	CellTensor U, dU_new, dU_old; 

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			U[i][j] = U_inlet;
		}
	}

	BoundaryConditions BoundaryTypes(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Symmetry);    

	RampGrid grid(Nx, Ny, 10, 10, 10, 6, 15); 

	string gridtype; 
	if (dynamic_cast<RampGrid*>(&grid)) gridtype = "Ramp";
	else if (dynamic_cast<CylinderGrid*>(&grid)) gridtype = "Cylinder";
	else gridtype = "Unknown"; 

	cout << "Running DPLR for " << Nx << " by " << Ny << " " << gridtype << "..." << endl;

	iFaceTensor i_Fluxes;
	jFaceTensor j_Fluxes;

	iFaceTesseract i_plus_Jacobians, i_minus_Jacobians;
	jFaceTesseract j_plus_Jacobians, j_minus_Jacobians; 

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < 4; ++k) {
				dU_new[i][j][k] = 0.0; 
				dU_old[i][j][k] = 0.0; 
			}
		}
	}

	double outer_residual = 1.0;

	while (outer_residual >= 0.000001) {

		dt = calculate_dt(U, grid, CFL);  

		solveOneTimestep(U, dU_new, U_inlet, dU_old, grid, BoundaryTypes, dt, CFL,
			i_Fluxes, j_Fluxes, i_plus_Jacobians, i_minus_Jacobians, j_plus_Jacobians, j_minus_Jacobians);


		outer_residual = calculateResidual(U, grid, i_Fluxes, j_Fluxes); 

		if (counter == 0) outer_residual = 1;		 
		counter++;   

		if (counter % 20 == 0) {
			auto end = TIME;
			DURATION duration = end - start; 
			cout << "Iteration: " << counter << "\t Residual: " << outer_residual << "\tElapsed time: " << duration.count() << endl;
		}

	}

	string filename = to_string(Nx) + "x" + to_string(Ny) + "_" + gridtype + "_Solution.csv"; 
	write2DCSV(filename, U, grid);   

	printMemoryUsage();

	return 0;
}