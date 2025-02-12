// 1D FVS Library.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "framework.h"
#include "1DFVSLibrary.h"

// Function that calculates dt based off of CFL
double calculate_dt(Matrix1D& U, const double& dx, const int& Nx, const double CFL) {

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

Matrix primtoCons(const Matrix& V) {
	double gamma = 1.4;

	Matrix U(3, Vector(1, 0.0));
	U[0][0] = V[0][0];
	U[1][0] = V[0][0] * V[1][0];
	U[2][0] = V[2][0] / (gamma - 1) + 0.5 * V[0][0] * V[1][0] * V[1][0];

	return U;
}

Matrix constoPrim(const Matrix& U) {
	double gamma = 1.4;

	Matrix V(3, Vector(1, 0.0));
	V[0][0] = U[0][0];
	V[1][0] = U[1][0] / U[0][0];
	V[2][0] = (U[2][0] - 0.5 * V[0][0] * V[1][0] * V[1][0]) * (gamma - 1);

	return V;
}

Matrix A_Minus(const Matrix& U, int nx) {
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

	Matrix Flux = dudv * K * Lambda * Kinv * dvdu;

	return Flux;
}

Matrix A_Plus(const Matrix& U, int nx) {
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

	Matrix Flux = dudv * K * Lambda * Kinv * dvdu;
	return Flux;
}

Duo Boundary1D(const Matrix U, int nx, const string type) {

	if (type == "inlet")
	{
		Matrix E = zeros(3, 3);
		return make_pair(U, E);
	}

	else if (type == "outlet") {
		Matrix E = ones(3, 3);
		return make_pair(U, E);
	}

	else {
		throw invalid_argument("Unknown boundary type: " + type);
	}


}

Matrix1D solveLine(Matrix1D U, const double dx, const int Nx, const double dt, Matrix A, Matrix B, Matrix C, Matrix F, Matrix alpha, Matrix2D v, Matrix2D g) {


	Matrix1D dU_new(3, Nx);

	// Define boundaries at the top and bottom of the line
	string leftBoundary = "outlet";
	string rightBoundary = "outlet";

	// Grab boundary values
	Duo leftBC = Boundary1D(U.vars(0), -1, leftBoundary);
	Duo rightBC = Boundary1D(U.vars(Nx - 1), 1, rightBoundary);

	// Top Boundary
	A = dx / dt * identity(3) + (A_Plus(U.vars(Nx - 1), 1) + rightBC.second * A_Minus(rightBC.first, 1)) + A_Plus(U.vars(Nx - 1), -1);
	C = A_Minus(U.vars(Nx - 2), -1);
	F = -1 * (A_Plus(U.vars(Nx - 1), 1) * U.vars(Nx - 1) + A_Minus(rightBC.first, 1) * rightBC.first + A_Plus(U.vars(Nx - 1), -1) * U.vars(Nx - 1) + A_Minus(U.vars(Nx - 2), -1) * U.vars(Nx - 2));

	alpha = A;
	v.strip(Nx - 1) = F / alpha;
	g.strip(Nx - 1) = C / alpha;

	// Middle Boundry
	for (int j = Nx - 2; j > 0; --j) {

		B = A_Minus(U.vars(j + 1), 1);
		A = dx / dt * identity(3) + A_Plus(U.vars(j), 1) + A_Plus(U.vars(j), -1);
		C = A_Minus(U.vars(j - 1), -1);
		F = -1 * (A_Plus(U.vars(j), 1) * U.vars(j) + A_Minus(U.vars(j + 1), 1) * U.vars(j + 1) + A_Plus(U.vars(j), -1) * U.vars(j) + A_Minus(U.vars(j - 1), -1) * U.vars(j - 1));

		alpha = A - B * g.strip(j + 1);
		g.strip(j) = C / alpha;
		v.strip(j) = (F - B * v.strip(j + 1)) / alpha;
	}

	//  Bottom boundary
	B = A_Minus(U.vars(1), 1);
	A = dx / dt * identity(3) + (A_Plus(U.vars(0), -1) + leftBC.second * A_Minus(leftBC.first, -1)) + A_Plus(U.vars(0), 1);
	F = -1 * (A_Plus(U.vars(0), 1) * U.vars(0) + A_Minus(U.vars(1), 1) * U.vars(1) + A_Plus(U.vars(0), -1) * U.vars(0) + A_Minus(leftBC.first, -1) * leftBC.first);


	alpha = A - B * g.strip(1);
	v.strip(0) = (F - B * v.strip(1)) / alpha;

	// Calculate dU
	dU_new.vars(0) = v.strip(0);

	for (int j = 1; j < Nx; ++j) {
		dU_new.vars(j) = v.strip(j) - g.strip(j) * dU_new.vars(j - 1);
	}

	return dU_new;
}

// Function that creates a .vtk file for Paraview postprocessing
void saveVectorsToCSV(const string& filename,
	const Vector& density,
	const Vector& u_vel,
	const Vector& pressure,
	const Vector& x) {

	ofstream csvFile(filename);

	if (!csvFile.is_open()) {
		cerr << "Error: Unable to open file " << filename << endl;
		return;
	}

	// Determine the number of points (assume all vectors have the same size)
	size_t numPoints = density.size();

	csvFile << "density, u velocity, pressure, x points" << endl;
	for (size_t i = 0; i < numPoints; ++i) {
		csvFile << density[i] << "," << u_vel[i] << "," << pressure[i] << "," << x[i] << endl;  // Use only x-component for simplicity
	}

	csvFile.close();
	cout << "VTK file saved to " << filename << endl;
}