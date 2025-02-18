#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <fstream> 
#include "GridGenerator.h" 
#include "2DFVSLibrary.h" 
#include "LinearAlgebra.h"
#include <omp.h> 

using namespace std;


#define TIME chrono::high_resolution_clock::now(); 
#define DURATION chrono::duration<double> duration;   

int main() {


	cout << "Running DPLR Finite Volume..." << endl;

	auto start = TIME;

	int Nx = 400, Ny = 200, n = 4, counter = 0;
	double CFL = 1.0, dt;

	double p = 10000.0, T = 300.0, R = 287.0, M = 2.5, a = sqrt(gamma * R * T), u = M * a, v = 0, rho = p / (R * T);

	Vector V_inlet = { rho, u, v, p };
	Vector U_inlet = primtoCons(V_inlet);

	Tensor U(Nx, Matrix(Ny, Vector(n, 0.0)));

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			U[i][j] = U_inlet;
		}
	}

	Tensor dU_new(Nx, Matrix(Ny, Vector(n, 0.0)));
	Tensor dU_old(Nx, Matrix(Ny, Vector(n, 0.0)));

	BoundaryConditions BoundaryTypes(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Symmetry);

	RampGrid grid(Nx, Ny, 10, 10, 10, 10, 15);

	Tensor i_Fluxes(Nx + 1, Matrix(Ny, Vector(n))), j_Fluxes(Nx, Matrix(Ny + 1, Vector(n)));

	Tesseract i_plus_Jacobians(Nx + 1, Tensor(Ny, Matrix(n, Vector(n)))), i_minus_Jacobians(Nx + 1, Tensor(Ny, Matrix(n, Vector(n)))),
		j_plus_Jacobians(Nx, Tensor(Ny + 1, Matrix(n, Vector(n)))), j_minus_Jacobians(Nx, Tensor(Ny + 1, Matrix(n, Vector(n))));

	double outer_residual = 10.0;

	while (outer_residual >= 1e-8) {


		solveOneTimestep(U, dU_new, U_inlet, dU_old, grid, BoundaryTypes, Nx, Ny, dt, CFL,
			i_Fluxes, j_Fluxes, i_plus_Jacobians, i_minus_Jacobians, j_plus_Jacobians, j_minus_Jacobians);

		outer_residual = calculateResidual(U, grid, Nx, Ny, i_Fluxes, j_Fluxes);
		if (counter == 0) outer_residual = 1;

		if (counter % 50 == 0) {
			auto end = TIME;
			DURATION duration = end - start;
			cout << "Iteration: " << counter << "\t Residual: " << outer_residual << "\t Elapsed time: " << duration.count() << endl;
		}
		counter++;
	}

	auto end = TIME;
	DURATION duration = end - start;
	cout << "Entire program took " << duration.count() << "seconds." << endl;

	writeVTK("2DRamp.csv", U, grid, Nx, Ny);
	writeCSV("1DRamp.csv", U, grid, 0, Nx);

	printMemoryUsage();

	return 0;
}