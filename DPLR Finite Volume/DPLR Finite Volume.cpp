// DPLR Finite Volume.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "LinearAlgebra.h"
#include "GridGenerator.h"
#include "2DFVSLibrary.h"
#include <chrono>


int main() {
	auto start1 = chrono::high_resolution_clock::now(); // Start the timer  
	cout << "DPLR Running" << endl << "--------------------------" << endl;

	int Nx = 100;
	int Ny = 50; 
	int Ramp_Angle = 15;

	RampGrid grid(Nx, Ny, 10, 10, 10, 7.5, Ramp_Angle);

	string filename = "Inviscid_Ramp_DPLR.vtk"; 

	BoundaryCondition leftBoundary = BoundaryCondition::Inlet;
	BoundaryCondition rightBoundary = BoundaryCondition::Outlet;
	BoundaryCondition bottomBoundary = BoundaryCondition::Symmetry;
	BoundaryCondition topBoundary = BoundaryCondition::Symmetry;

	double M = 2.5;
	double P = 10000;
	double T = 300;

	double R = 287;
	double a = sqrt(gamma * R * T);
	double rho = P / (R * T); 
	double u = M * a;
	double v = 0;
	double CFL = 1.0;

	Matrix V_inlet = { {rho}, {u}, {v}, {P} };
	Matrix U_inlet = primtoCons(V_inlet);

	Matrix2D U(4, Nx, Ny), U_old(4, Nx, Ny), dU(4, Nx, Ny), dU_new(4, Nx, Ny), dU_old(4, Nx, Ny);

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			U.vars(i, j) = U_inlet; 
		}
	}

	double outer_residual = 1.0;
	int n = 1;
	double t_old = 0.0;
	Vector t;
	double t0 = 0.0;
	t.push_back(t0);  

	double dt; 

	while (outer_residual >= 1e-6) {

		dt = calculate_dt(U, grid, Nx, Ny, CFL); 
		dU_old = dU;
		U_old = U;

		dU = solveOneTimestep(U, dU_old, U_inlet, grid, dt, Nx, Ny, leftBoundary, rightBoundary, bottomBoundary, topBoundary);
		dU.NanorInf(); 

		U += dU; 

		t.push_back(t.back() + dt); 
		outer_residual = calculateResidual(U, grid, Nx, Ny);
		
		if (n % 20 == 0) {
			auto end1 = chrono::high_resolution_clock::now(); // Stop the timer 
			chrono::duration<double> duration1 = end1 - start1; // Calculate the duration 
			cout << "Iteration " << n << ":\t Residual = " << outer_residual << "\t dt = " << dt << "\t Elapsed time = " << duration1.count() << endl;  
			writeVTK(filename, U, grid, Nx, Ny); 
		}		
		n++; 		
	}

	writeVTK(filename, U, grid, Nx, Ny); 

	cout << endl;
	auto end1 = chrono::high_resolution_clock::now(); // Stop the timer 
	chrono::duration<double> duration1 = end1 - start1; // Calculate the duration 
	cout << "DPLR took " << duration1.count() << " seconds to complete." << endl;

	return 0;


}

