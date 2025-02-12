// 1D Sod Shock Tube (With Libraries).cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "LinearAlgebra.h"
#include "1DFVSLibrary.h"


using namespace std;
constexpr double CFL = 0.3;

double calculate_dt(Matrix1D& U, const double& dx, const int& Nx) {

	double dt = 0, dtold = 0;
	Matrix V;
	double a;
	for (int i = 0; i < Nx; ++i) {

		V = constoPrim(U.vars(i));
		a = sqrt(gamma * V[2][0] / V[0][0]);


		dt = CFL * dx / (fabs(V[1][0]) + a);

		if (dt > dtold) {
			dtold = dt;
		}

	}

	return dtold;
}

int main() {

	auto start1 = chrono::high_resolution_clock::now(); // Start the timer 

	cout << "1D Sod Shock with library (.vars method)" << endl;

	// Spatial variables
	double rhol = 1.0, pl = 1.0, ul = 0.0, rhor = 0.125, pr = 0.1, ur = 0.0;
	const int Nx = 1000;
	const double Lx = 2;
	const double dx = Lx / (Nx - 1);

	// Temporal variables 
	double t = 0.0, tfinal = 0.5;
	double dt = CFL * dx;
	int n = 0;

	// x values for plotting
	Vector x(Nx, 0.0);

	for (int i = 0; i < Nx; ++i) {
		x[i] = -1 + i * dx;
	}


	// Initial primitive variables
	Matrix VL = { {rhol}, {ul}, {pl} };
	Matrix VR = { {rhor}, {ur}, {pr} };

	// Initial conserved variables. A function is called to change primitives to conserved.
	Matrix UL = primtoCons(VL);
	Matrix UR = primtoCons(VR);

	// Initialize empty U matrix
	Matrix1D U(3, Nx);
	Matrix1D Unew(3, Nx);

	// Fill in U for t = 0
	for (int i = 0; i < Nx; ++i) {
		if (i <= Nx / 2) {
			U.vars(i) = UL;
		}
		else {
			U.vars(i) = UR;
		}
	}

	Matrix FR, FL;

	// Face normals for 1D
	int nxn = -1;
	int nxp = 1;



	// Main while loop 
	while (t <= tfinal) {

		// Loop through all cells

		dt = calculate_dt(U, dx, Nx);
		for (int i = 1; i < Nx - 1; ++i) {


			// All fluxes
			FR = A_Plus(U.vars(i), nxp) * U.vars(i) + A_Minus(U.vars(i + 1), nxp) * U.vars(i + 1);
			FL = A_Plus(U.vars(i), nxn) * U.vars(i) + A_Minus(U.vars(i - 1), nxn) * U.vars(i - 1);

			// Update U and check if NaN or inf
			Unew.vars(i) = U.vars(i) - (dt / dx * (FR + FL));  

			//for (int k = 0; k < 3; ++k) {
			//	Unew[k][i] = U[k][i] - dt / dx * (FR[k][0] + FL[k][0]); 
			//}
		}

		U = Unew; 

		// Make sure BC are enforced		
		U.vars(0) = UL;
		U.vars(Nx - 1) = UR;

		// Update t and n counter
		t += dt;
		n++;

		// Look at progress

		if (n % 100 == 0) {
			cout << "n = " << n << "\t" << "t = " << t << endl;
		}
	}



	cout << endl;
	auto end1 = chrono::high_resolution_clock::now(); // Stop the timer 
	chrono::duration<double> duration1 = end1 - start1; // Calculate the duration 
	cout << "Sod shock tube WITH library took " << duration1.count() << " seconds to complete." << endl;

	// Create vectors of individual variables for plotting
	Vector density(Nx), u_vel(Nx), pressure(Nx);
	Matrix Primitives;

	for (int i = 0; i < Nx; ++i) {

		Primitives = constoPrim(U.vars(i));
		density[i] = Primitives[0][0];
		u_vel[i] = Primitives[1][0];
		pressure[i] = Primitives[2][0];
	}

	saveVectorsToCSV("1DSodShock.csv", density, u_vel, pressure, x);

	cout << "t = " << t << endl;

	return 0;
}