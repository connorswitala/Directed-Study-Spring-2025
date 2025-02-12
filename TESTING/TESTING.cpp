
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <fstream> 
#include <chrono>
#include <omp.h>

using namespace std;


typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

// Declar functions

Matrix primtoCons(Matrix& V);
Matrix constoPrim(Matrix& U);
Matrix Flux_Plus(Matrix& U, int nx);
Matrix Flux_Minus(Matrix& U, int nx);
void saveVectorsToCSV(const std::string& filename,
	const std::vector<double>& vector1,
	const std::vector<double>& vector2,
	const std::vector<double>& vector3,
	const std::vector<double>& vector4);




int main() {
	omp_set_num_threads(16);
	auto start1 = chrono::high_resolution_clock::now(); // Start the timer  

	// Spatial variables
	double rhol = 1.0, pl = 1.0, ul = 0.0, rhor = 0.125, pr = 0.1, ur = 0.0, gamma = 1.4;
	const int Nx = 500;
	const double Lx = 2;
	const double dx = Lx / (Nx - 1);

	// Temporal variables 
	double t = 0.0, tfinal = 0.3, CFL = 0.1;
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
	Matrix U(3, Vector(Nx, 0.0));
	Matrix U_new(3, Vector(Nx, 0.0));

	// Fill in U for t = 0
	for (int i = 0; i < Nx; ++i) {
		if (i <= Nx / 2) {
			U[0][i] = UL[0][0];
			U[1][i] = UL[1][0];
			U[2][i] = UL[2][0];
		}
		else {
			U[0][i] = UR[0][0];
			U[1][i] = UR[1][0];
			U[2][i] = UR[2][0];
		}
	}

	// Face normals for 1D
	int nxn = -1;
	int nxp = 1;

	// Main while loop 
	while (t <= tfinal) {
	#pragma omp parallel
			{
				// Declare temporary variables ONCE per thread
				Matrix U_i_minus_one, U_i, U_i_plus_one;
				Matrix FPR, FPL, FMR, FML; 

				// Loop through all cells
				#pragma omp parallel for 
				for (int i = 1; i < Nx - 1; ++i) {
					// Create temporary matrix of U to pass into functions and get flux
					U_i_minus_one = { {U[0][i - 1]}, {U[1][i - 1]}, {U[2][i - 1]} };
					U_i = { {U[0][i]}, {U[1][i]}, {U[2][i]} };
					U_i_plus_one = { {U[0][i + 1]}, {U[1][i + 1]}, {U[2][i + 1]} };

					// All fluxes
					FPR = Flux_Plus(U_i, nxp);
					FPL = Flux_Plus(U_i, nxn);
					FMR = Flux_Minus(U_i_plus_one, nxp);
					FML = Flux_Minus(U_i_minus_one, nxn);

					// Update U_new and check if NaN or inf
					for (int j = 0; j < 3; ++j) {
						U_new[j][i] = U[j][i] - dt / dx * (FPR[j][0] + FMR[j][0] + FPL[j][0] + FML[j][0]);
					}

				}
			}

		for (int i = 1; i < Nx - 1; ++i) {
			for (int k = 0; k < 3; ++k) {
				U[k][i] = U_new[k][i];
			}
 		}

		// Make sure BC are enforced		
		U[0][0] = UL[0][0];
		U[1][0] = UL[1][0];
		U[2][0] = UL[2][0];

		U[0][Nx - 1] = UR[0][0];
		U[1][Nx - 1] = UR[1][0];
		U[2][Nx - 1] = UR[2][0];


		// Update t and n counter
		t += dt;
		n++;

		// Look at progress

		if (n % 100 == 0) {
			cout << "n = " << n << "\t" << "t = " << t << endl;
		}
	}

	auto end1 = chrono::high_resolution_clock::now(); // Stop the timer
	chrono::duration<double> duration1 = end1 - start1; // Calculate the duration
	cout << "Sod shock tube WITH library took " << duration1.count() << " seconds to complete." << endl;

	// Create vectors of individual variables for plotting
	Vector density(Nx), u_vel(Nx), pressure(Nx);
	Matrix U_pass, Primitives;

	for (int i = 0; i < Nx; ++i) {
		U_pass = { {U[0][i]}, {U[1][i]}, {U[2][i]} };

		Primitives = constoPrim(U_pass);
		density[i] = Primitives[0][0];
		u_vel[i] = Primitives[1][0];
		pressure[i] = Primitives[2][0];
	}

	saveVectorsToCSV("1DSodShock.csv", density, u_vel, pressure, x);

	cout << "t = " << t << endl;


	return 0;
}

// Function that multiplies matrices
Matrix MultiplyMatrices(Matrix& A, Matrix& B) {

	size_t Arows = A.size();
	size_t Brows = B.size();
	size_t Acols = A[0].size();
	size_t Bcols = B[0].size();

	if (Acols != Brows) {
		cout << "A is a " << Arows << " by " << Acols << " matrix." << endl;
		cout << "B is a " << Brows << " by " << Bcols << " matrix." << endl;
		cout << "Improper dimensions for matrix multiplication!" << endl << endl;
		exit(1);
	}

	Matrix Resultant_Matrix(Arows, vector<double>(Bcols, 0.0));

	for (size_t i = 0; i < Arows; ++i) {
		for (size_t j = 0; j < Bcols; ++j) {
			for (size_t k = 0; k < Acols; ++k) {
				Resultant_Matrix[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	return Resultant_Matrix;
}
// Function that converts primitives to conserved variables
Matrix primtoCons(Matrix& V) {
	double gamma = 1.4;

	Matrix U(3, vector<double>(1, 0.0));
	U[0][0] = V[0][0];
	U[1][0] = V[0][0] * V[1][0];
	U[2][0] = V[2][0] / (gamma - 1) + 0.5 * V[0][0] * V[1][0] * V[1][0];

	return U;
}
// Function that converts conserved variables to primitives
Matrix constoPrim(Matrix& U) {
	double gamma = 1.4;

	Matrix V(3, vector<double>(1, 0.0));
	V[0][0] = U[0][0];
	V[1][0] = U[1][0] / U[0][0];
	V[2][0] = (U[2][0] - 0.5 * V[0][0] * V[1][0] * V[1][0]) * (gamma - 1);

	return V;
}
// Function that calculates F+
Matrix Flux_Plus(Matrix& U, int nx) {
	double gamma = 1.4;

	Matrix V = constoPrim(U);
	double rho = V[0][0];
	double u = V[1][0];
	double p = V[2][0];

	double a = sqrt(gamma * p / rho);

	Matrix Lambda = {
		{0.5 * ((u * nx - a) + fabs(u * nx - a)), 0, 0 },
		{0, 0.5 * (u * nx + fabs(u * nx)), 0},
		{0, 0, 0.5 * ((u * nx + a) + fabs(u * nx + a)) }
	};

	Matrix K{
		{rho / (gamma * p), 1, rho / (gamma * p)},
		{-nx / sqrt(rho * gamma * p), 0, nx / sqrt(rho * gamma * p)},
		{1, 0, 1}
	};

	Matrix Kinv = {
		{0, -0.5 * nx * sqrt(gamma * p * rho), 0.5},
		{1, 0, -rho / (gamma * p)},
		{0, 0.5 * nx * sqrt(gamma * rho * p), 0.5}
	};

	Matrix dvdu = {
		{1, 0, 0},
		{-u / rho, 1 / rho, 0},
		{0.5 * u * u * (gamma - 1), -u * (gamma - 1), gamma - 1 }
	};

	Matrix dudv = {
		{1, 0, 0},
		{u, rho, 0},
		{0.5 * u * u, rho * u, 1 / (gamma - 1) }
	};

	Matrix int1 = MultiplyMatrices(dudv, K);
	Matrix int2 = MultiplyMatrices(int1, Lambda);
	Matrix int3 = MultiplyMatrices(int2, Kinv);
	Matrix int4 = MultiplyMatrices(int3, dvdu);
	Matrix Flux = MultiplyMatrices(int4, U);

	return Flux;
}
// Function that calculates F-
Matrix Flux_Minus(Matrix& U, int nx) {
	double gamma = 1.4;

	Matrix V = constoPrim(U);
	double rho = V[0][0];
	double u = V[1][0];
	double p = V[2][0];

	double a = sqrt(gamma * p / rho);

	Matrix Lambda = {
		{0.5 * ((u * nx - a) - fabs(u * nx - a)), 0, 0 },
		{0, 0.5 * (u * nx - fabs(u * nx)), 0},
		{0, 0, 0.5 * ((u * nx + a) - fabs(u * nx + a)) }
	};

	Matrix K{
		{rho / (gamma * p), 1, rho / (gamma * p)},
		{-nx / sqrt(rho * gamma * p), 0, nx / sqrt(rho * gamma * p)},
		{1, 0, 1}
	};

	Matrix Kinv = {
		{0, -0.5 * nx * sqrt(gamma * p * rho), 0.5},
		{1, 0, -rho / (gamma * p)},
		{0, 0.5 * nx * sqrt(gamma * rho * p), 0.5}
	};

	Matrix dvdu = {
		{1, 0, 0},
		{-u / rho, 1 / rho, 0},
		{0.5 * u * u * (gamma - 1), -u * (gamma - 1), gamma - 1 }
	};

	Matrix dudv = {
		{1, 0, 0},
		{u, rho, 0},
		{0.5 * u * u, rho * u, 1 / (gamma - 1) }
	};

	Matrix int1 = MultiplyMatrices(dudv, K);
	Matrix int2 = MultiplyMatrices(int1, Lambda);
	Matrix int3 = MultiplyMatrices(int2, Kinv);
	Matrix int4 = MultiplyMatrices(int3, dvdu);
	Matrix Flux = MultiplyMatrices(int4, U);



	return Flux;
}
// Function that creates a .vtk file for Paraview postprocessing
void saveVectorsToCSV(const string& filename,
	const std::vector<double>& density,
	const std::vector<double>& u_vel,
	const std::vector<double>& pressure,
	const std::vector<double>& x) {
	std::ofstream csvFile(filename);

	if (!csvFile.is_open()) {
		std::cerr << "Error: Unable to open file " << filename << std::endl;
		return;
	}


	// Determine the number of points (assume all vectors have the same size)
	size_t numPoints = density.size();

	csvFile << "density, u velocity, pressure, x points" << endl;
	for (size_t i = 0; i < numPoints; ++i) {
		csvFile << density[i] << "," << u_vel[i] << "," << pressure[i] << "," << x[i] << endl;  // Use only x-component for simplicity
	}

	csvFile.close();
	std::cout << "VTK file saved to " << filename << std::endl;
}