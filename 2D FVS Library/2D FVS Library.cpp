// 2D FVS Library.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "framework.h"
#include "2DFVSLibrary.h"
#include <omp.h>

// Function that converts primitives to conserved variables
Matrix primtoCons(const Matrix& V) {

	Matrix U(4, vector<double>(1, 0.0));
	U[0][0] = V[0][0];
	U[1][0] = V[0][0] * V[1][0];
	U[2][0] = V[0][0] * V[2][0];
	U[3][0] = V[3][0] / (gamma - 1) + 0.5 * V[0][0] * (V[1][0] * V[1][0] + V[2][0] * V[2][0]);

	return U;
}

// Function that converts conserved variables to primitives
Matrix constoPrim(const Matrix& U) {

	Matrix V(4, vector<double>(1, 0.0));
	V[0][0] = U[0][0];
	V[1][0] = U[1][0] / U[0][0];
	V[2][0] = U[2][0] / U[0][0];
	V[3][0] = (U[3][0] - 0.5 * V[0][0] * (V[1][0] * V[1][0] + V[2][0] * V[2][0])) * (gamma - 1);

	return V;
}

// Function that calculates ghost cells
Duo Boundary2D(BoundaryCondition type, const Matrix U, const Matrix U_inlet, const Point normals) {

	Matrix E, Ug, inside_Primitives; 
	Matrix ghost_Primitives(4, Vector(1, 0.0)); 
	double u = 0.0, v = 0.0; 
	switch (type) {

	case BoundaryCondition::Inlet: 
		E = zeros(4, 4); 
		return make_pair(U_inlet, E);  
	case BoundaryCondition::Outlet: 
		E = ones(4, 4); 
		return make_pair(U, E); 
	case BoundaryCondition::Wall: 
		Ug = ones(4, 1); 
		E = ones(4, 4); 
		
		return make_pair(Ug, E); 

	case BoundaryCondition::Symmetry: 

		inside_Primitives = constoPrim(U); 

		u = inside_Primitives[1][0];
		v = inside_Primitives[2][0];

		ghost_Primitives[0][0] = inside_Primitives[0][0]; 
		ghost_Primitives[1][0] = u - 2 * (u * normals.x + v * normals.y) * normals.x; 
		ghost_Primitives[2][0] = v - 2 * (u * normals.x + v * normals.y) * normals.y; 
		ghost_Primitives[3][0] = inside_Primitives[3][0]; 

		Ug = primtoCons(ghost_Primitives);

		E = { {1, 0, 0, 0},
					 {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
					 { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
					 { 0, 0, 0, 1 }
		};

		return make_pair(Ug, E);

	default:
		throw invalid_argument("Unknown boundary condition type.");
	}

}


// Writes to VTK file (Contout plot)
void writeVTK(const string& filename,
	Matrix2D& U, 
	const Grid& grid,
	const int Nx,
	const int Ny) {

	Matrix density(Nx, vector<double>(Ny, 0.0)), u_vel(Nx, vector<double>(Ny, 0.0)), v_vel(Nx, vector<double>(Ny, 0.0)), pressure(Nx, vector<double>(Ny, 0.0)); 
	Matrix Primitives(4, Vector(1, 0.0));

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) { 

			Primitives = constoPrim(U.vars(i, j)); 
			density[i][j] = Primitives[0][0];
			u_vel[i][j] = Primitives[1][0];
			v_vel[i][j] = Primitives[2][0];
			pressure[i][j] = Primitives[3][0];
		}
	}
	ofstream file(filename);

	// Header
	file << "# vtk DataFile Version 3.0\n";
	file << "2D CFD data\n";
	file << "ASCII\n";
	file << "DATASET STRUCTURED_GRID\n";
	file << "DIMENSIONS " << Nx << " " << Ny << " 1\n";	// Add assertions to check vector sizes


	// Write points
	file << "POINTS " << Nx * Ny << " float\n";
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			file << grid.Center(i, j).x << " " << grid.Center(i, j).y << " 0.0\n";
		}
	}


	// Point data (scalars and vectors)
	file << "\nPOINT_DATA " << Nx * Ny << "\n";

	// Density
	file << "SCALARS density float\n";
	file << "LOOKUP_TABLE default\n";
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			file << density[i][j] << "\n";
		}
	}


	// Pressure
	file << "\nSCALARS pressure float\n";
	file << "LOOKUP_TABLE default\n";
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			file << pressure[i][j] << "\n";
		}
	}

	// Velocity (as vectors)
	file << "\nVECTORS velocity float\n";
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			file << u_vel[i][j] << " " << v_vel[i][j] << " 0.0\n"; 
		}
	}

	file.close();
}

// Writes to CSV file (Line plot)
void writeCSV(const string& filename,
	const Matrix& density,
	const Matrix& u_velocity,
	const Matrix& v_velocity,
	const Matrix& pressure,
	const Grid& grid,
	const int n,
	const int Nx) {

	ofstream file(filename);

	file << "density, u_velocity, v_velocity, pressure, x_points, y_points" << endl;

	for (int i = 0; i < Nx; ++i) {
		file << density[i][n] << ", " << u_velocity[i][n] << ", " << v_velocity[i][n] << ", " << pressure[i][n] << ", " << grid.Center(i, n).x << ", " << grid.Center(i, n).y << ", " << "0.0\n" << endl;
	}

	file.close();

}


// Function that calculates A-
Matrix A_Minus(const Matrix& U, Point normals) {

	double nx = normals.x;
	double ny = normals.y;

	Matrix V = constoPrim(U);
	double rho = V[0][0];
	double u = V[1][0];
	double v = V[2][0];
	double p = V[3][0];

	double a = sqrt(gamma * p / rho);

	double uprime = u * nx + v * ny;
	double l1 = 0.5 * (uprime - a - fabs(uprime - a));
	double l2 = 0.5 * (uprime - fabs(uprime));
	double l3 = 0.5 * (uprime - fabs(uprime));
	double l4 = 0.5 * (uprime + a - fabs(uprime + a));


	Matrix M = {
		{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)},
		{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
	};

	Matrix dvdu = {
		{1, 0, 0, 0},
		{-u / rho, 1 / rho, 0, 0},
		{-v / rho, 0, 1 / rho, 0},
		{0.5 * (gamma - 1) * (u * u + v * v), -u * (gamma - 1), -v * (gamma - 1), gamma - 1}
	};

	Matrix dudv = {
		{1, 0, 0, 0},
		{u, rho, 0, 0},
		{v, 0, rho, 0},
		{0.5 * (u * u + v * v), rho * u, rho * v, 1 / (gamma - 1)}
	};

	Matrix A_Minus = dudv * M * dvdu;
	return A_Minus;

}

// Function that calculates A+
Matrix A_Plus(const Matrix& U, Point normals) {

	double nx = normals.x;
	double ny = normals.y;

	Matrix V = constoPrim(U);
	double rho = V[0][0];
	double u = V[1][0];
	double v = V[2][0];
	double p = V[3][0];

	double a = sqrt(gamma * p / rho);

	double uprime = u * nx + v * ny;
	double l1 = 0.5 * (uprime - a + fabs(uprime - a));
	double l2 = 0.5 * (uprime + fabs(uprime));
	double l3 = 0.5 * (uprime + fabs(uprime));
	double l4 = 0.5 * (uprime + a + fabs(uprime + a));


	Matrix M = {
		{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)},
		{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
	};

	Matrix dvdu = {
		{1, 0, 0, 0},
		{-u / rho, 1 / rho, 0, 0},
		{-v / rho, 0, 1 / rho, 0},
		{0.5 * (gamma - 1) * (u * u + v * v), -u * (gamma - 1), -v * (gamma - 1), gamma - 1}
	};

	Matrix dudv = {
		{1, 0, 0, 0},
		{u, rho, 0, 0},
		{v, 0, rho, 0},
		{0.5 * (u * u + v * v), rho * u, rho * v, 1 / (gamma - 1)}
	};

	Matrix A_Plus = dudv * M * dvdu;
	return A_Plus;
}

// Function that calculates dt
double calculate_dt(Matrix2D& U, const Grid& grid, int Nx, int Ny, double CFL) {

	Matrix V(4, Vector(1, 0.0));
	double dx, dy, c, dt_old, dt = 1;

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			V = constoPrim(U.vars(i, j));
			dx = min(grid.BArea(i, j), grid.TArea(i, j));
			dy = min(grid.LArea(i, j), grid.RArea(i, j));
			c = sqrt(gamma * V[3][0] / V[0][0]);

			dt_old = CFL / (fabs(V[1][0] / dx) + fabs(V[2][0] / dy) + c * sqrt(1 / (dx * dx) + 1 / (dy * dy)));

			if (dt_old < dt) dt = dt_old;

		}
	}
	return dt;
}

double calculateInnerResidual(Matrix2D& U, Matrix2D& dU_old, Matrix2D& dU_new, const Grid& grid, const int Nx, const int Ny, const double dt, int i) {

	double lineresidual = 0.0;
	Matrix res(4, Vector(1, 0.0));

	Matrix A; // Relates to U(j+1) 
	Matrix B; // Relates to U(j) 
	Matrix C; // Relates to U(j-1) 
	Matrix F; // Right hand side 


	// Middle Boundary
	for (int j = Ny - 2; j > 0; --j) {

		B = A_Minus(U.vars(i, j + 1), grid.TNorms(i, j)) * grid.TArea(i, j); 

		A = grid.Volume(i, j) / dt * identity(4)
			+ A_Plus(U.vars(i, j), grid.LNorms(i, j)) * grid.LArea(i, j) 
			+ A_Plus(U.vars(i, j), grid.RNorms(i, j)) * grid.RArea(i, j)
			+ A_Plus(U.vars(i, j), grid.BNorms(i, j)) * grid.BArea(i, j)
			+ A_Plus(U.vars(i, j), grid.TNorms(i, j)) * grid.TArea(i, j);

		C = A_Minus(U.vars(i, j - 1), grid.BNorms(i, j)) * grid.BArea(i, j);

		F = -1 * ((A_Plus(U.vars(i, j), grid.LNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i - 1, j), grid.LNorms(i, j)) * U.vars(i - 1, j)) * grid.LArea(i, j)			// Left Side		
			+ (A_Plus(U.vars(i, j), grid.RNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i + 1, j), grid.RNorms(i, j)) * U.vars(i + 1, j)) * grid.RArea(i, j)		// Right Side
			+ (A_Plus(U.vars(i, j), grid.BNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i, j - 1), grid.BNorms(i, j)) * U.vars(i, j - 1)) * grid.BArea(i, j)			// Bottom Side
			+ (A_Plus(U.vars(i, j), grid.TNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i, j + 1), grid.TNorms(i, j)) * U.vars(i, j + 1)) * grid.TArea(i, j)  // Top boundary 
			+ A_Minus(U.vars(i - 1, j), grid.LNorms(i, j)) * dU_old.vars(i - 1, j) * grid.LArea(i, j)													// Left dU addition
			+ A_Minus(U.vars(i + 1, j), grid.RNorms(i, j)) * dU_old.vars(i + 1, j) * grid.RArea(i, j)													// Right dU addition 
			);

		res = B * dU_new.vars(i, j + 1) + A * dU_new.vars(i, j) + C * dU_new.vars(i, j - 1) - F;
		lineresidual = lineresidual + res[0][0] * res[0][0];
	}

	return lineresidual; 

}

double calculateResidual(Matrix2D& U, const Grid& grid, int Nx, int Ny) {

	double res = 0.0; 
	Matrix intres(4, Vector(1, 0.0)); 

	for (int i = 1; i < Nx - 1; ++i) {
		for (int j = 1; j < Ny - 1; ++j) {
			intres = (A_Plus(U.vars(i, j), grid.LNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i - 1, j), grid.LNorms(i, j)) * U.vars(i - 1, j)) * grid.LArea(i, j)			// Left Side		
				+ (A_Plus(U.vars(i, j), grid.RNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i + 1, j), grid.RNorms(i, j)) * U.vars(i + 1, j)) * grid.RArea(i, j)		// Right Side
				+ (A_Plus(U.vars(i, j), grid.BNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i, j - 1), grid.BNorms(i, j)) * U.vars(i, j - 1)) * grid.BArea(i, j)			// Bottom Side
				+ (A_Plus(U.vars(i, j), grid.TNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i, j + 1), grid.TNorms(i, j)) * U.vars(i, j + 1)) * grid.TArea(i, j);   

			res = res + intres[0][0] / grid.Volume(i,j) * intres[0][0] / grid.Volume(i, j);
		}
	}

	res = sqrt(res); 
	return res; 
	

}

Matrix2D solveOneTimestep(Matrix2D& U, Matrix2D& dU_old, Matrix U_inlet, const Grid& grid, double dt, int Nx, int Ny, BoundaryCondition left, BoundaryCondition right, BoundaryCondition bottom, BoundaryCondition top) {

	double inner_residual = 1.0;
	Matrix2D dU_new(4, Nx, Ny); 

	while (inner_residual >= 1e-8) {
		
		inner_residual = 0.0;
			
		dU_new.line(0) = solveLeftLine(U, U_inlet, dU_old, grid, left, top, bottom, 0, Nx, Ny, dt);  

		
		for (int i = 1; i < Nx - 1; ++i) {
			dU_new.line(i) = solveMiddleLine(U, U_inlet, dU_old, grid, top, bottom, i, Nx, Ny, dt);
		}	
	
		dU_new.line(Nx - 1) = solveRightLine(U, U_inlet, dU_old, grid, right, top, bottom, Nx - 1, Nx, Ny, dt);    	

		for (int i = 1; i < Nx - 1; ++i) {
			inner_residual = inner_residual + calculateInnerResidual(U, dU_old, dU_new, grid, Nx, Ny, dt, i);
		} 

		inner_residual = sqrt(inner_residual);
		swap(dU_old, dU_new);  
	}

	return dU_old; 

}


// Solve middle lines using DPLR
Matrix1D solveMiddleLine(Matrix2D& U, Matrix U_inlet, Matrix2D& dU_old, const Grid& grid, 
	BoundaryCondition toptype, BoundaryCondition bottomtype, 
	const int i, const int Nx, const int Ny, const double dt) { 

	Matrix A; // Relates to U(j+1) 
	Matrix B; // Relates to U(j) 
	Matrix C; // Relates to U(j-1) 
	Matrix F; // Right hand side 

	// Intermediate matrices
	Matrix alpha;
	Matrix2D v(4, 1, Nx);
	Matrix2D g(4, 4, Nx);

	Matrix1D dU(4, Ny);

	// Grab boundary values
	Duo bottomBC = Boundary2D(bottomtype, U.vars(i, 0), U_inlet, grid.BNorms(i, 0)); 
	Duo topBC = Boundary2D(toptype, U.vars(i, Ny - 1), U_inlet, grid.TNorms(i, Ny - 1));    


	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity(4)
		+ (A_Plus(U.vars(i, Ny - 1), grid.TNorms(i, Ny - 1)) + topBC.second * A_Minus(topBC.first, grid.TNorms(i, Ny - 1))) * grid.TArea(i, Ny - 1)
		+ A_Plus(U.vars(i, Ny - 1), grid.RNorms(i, Ny - 1)) * grid.RArea(i, Ny - 1)
		+ A_Plus(U.vars(i, Ny - 1), grid.BNorms(i, Ny - 1)) * grid.BArea(i, Ny - 1)
		+ A_Plus(U.vars(i, Ny - 1), grid.LNorms(i, Ny - 1)) * grid.LArea(i, Ny - 1); 

	C = A_Minus(U.vars(i, Ny - 2), grid.BNorms(i, Ny - 1)) * grid.BArea(i, Ny - 1);

	F = -1 * ((A_Plus(U.vars(i, Ny - 1), grid.LNorms(i, Ny - 1)) * U.vars(i, Ny - 1) + A_Minus(U.vars(i - 1, Ny - 1), grid.LNorms(i, Ny - 1)) * U.vars(i - 1, Ny - 1)) * grid.LArea(i, Ny - 1)	// Left Side		
		+ (A_Plus(U.vars(i, Ny - 1), grid.RNorms(i, Ny - 1)) * U.vars(i, Ny - 1) + A_Minus(U.vars(i + 1, Ny - 1), grid.RNorms(i, Ny - 1)) * U.vars(i + 1, Ny - 1)) * grid.RArea(i, Ny - 1)	// Right Side
		+ (A_Plus(U.vars(i, Ny - 1), grid.BNorms(i, Ny - 1)) * U.vars(i, Ny - 1) + A_Minus(U.vars(i, Ny - 2), grid.BNorms(i, Ny - 1)) * U.vars(i, Ny - 2)) * grid.BArea(i, Ny - 1)			// Bottom Side
		+ (A_Plus(U.vars(i, Ny - 1), grid.TNorms(i, Ny - 1)) * U.vars(i, Ny - 1) + A_Minus(topBC.first, grid.TNorms(i, Ny - 1)) * topBC.first) * grid.TArea(i, Ny - 1)						// Top boundary
		+ A_Minus(U.vars(i - 1, Ny - 1), grid.LNorms(i, Ny - 1)) * dU_old.vars(i - 1, Ny - 1) * grid.LArea(i, Ny - 1)																		// Left dU addition
		+ A_Minus(U.vars(i + 1, Ny - 1), grid.RNorms(i, Ny - 1)) * dU_old.vars(i + 1, Ny - 1) * grid.RArea(i, Ny - 1)																		// Right dU addition
		);

	alpha = A;
	v.strip(Ny - 1) = F / alpha;
	g.strip(Ny - 1) = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		B = A_Minus(U.vars(i, j + 1), grid.TNorms(i, j)) * grid.TArea(i, j);

		A = grid.Volume(i, j) / dt * identity(4)
			+ A_Plus(U.vars(i, j), grid.LNorms(i, j)) * grid.LArea(i, j)
			+ A_Plus(U.vars(i, j), grid.RNorms(i, j)) * grid.RArea(i, j)
			+ A_Plus(U.vars(i, j), grid.BNorms(i, j)) * grid.BArea(i, j)
			+ A_Plus(U.vars(i, j), grid.TNorms(i, j)) * grid.TArea(i, j);

		C = A_Minus(U.vars(i, j - 1), grid.BNorms(i, j)) * grid.BArea(i, j);

		F = -1 * ((A_Plus(U.vars(i, j), grid.LNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i - 1, j), grid.LNorms(i, j)) * U.vars(i - 1, j)) * grid.LArea(i, j)			// Left Side		
			+ (A_Plus(U.vars(i, j), grid.RNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i + 1, j), grid.RNorms(i, j)) * U.vars(i + 1, j)) * grid.RArea(i, j)		// Right Side
			+ (A_Plus(U.vars(i, j), grid.BNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i, j - 1), grid.BNorms(i, j)) * U.vars(i, j - 1)) * grid.BArea(i, j)			// Bottom Side
			+ (A_Plus(U.vars(i, j), grid.TNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i, j + 1), grid.TNorms(i, j)) * U.vars(i, j + 1)) * grid.TArea(i, j)  // Top boundary 
			+ A_Minus(U.vars(i - 1, j), grid.LNorms(i, j)) * dU_old.vars(i - 1, j) * grid.LArea(i, j)													// Left dU addition
			+ A_Minus(U.vars(i + 1, j), grid.RNorms(i, j)) * dU_old.vars(i + 1, j) * grid.RArea(i, j)													// Right dU addition
			);

		alpha = A - B * g.strip(j + 1);
		g.strip(j) = C / alpha;
		v.strip(j) = (F - B * v.strip(j + 1)) / alpha;
	}

	//  Bottom boundary
	B = A_Minus(U.vars(i, 1), grid.TNorms(i, 0)) * grid.TArea(i, 0);

	A = grid.Volume(i, 0) / dt * identity(4)
		+ A_Plus(U.vars(i, 0), grid.LNorms(i, 0)) * grid.LArea(i, 0)
		+ A_Plus(U.vars(i, 0), grid.RNorms(i, 0)) * grid.RArea(i, 0)
		+ A_Plus(U.vars(i, 0), grid.BNorms(i, 0)) * grid.BArea(i, 0)
		+ A_Plus(U.vars(i, 0), grid.TNorms(i, 0)) * grid.TArea(i, 0);

	F = -1 * ((A_Plus(U.vars(i, 0), grid.LNorms(i, 0)) * U.vars(i, 0) + A_Minus(U.vars(i - 1, 0), grid.LNorms(i, 0)) * U.vars(i - 1, 0)) * grid.LArea(i, 0)		// Left Side		
		+ (A_Plus(U.vars(i, 0), grid.RNorms(i, 0)) * U.vars(i, 0) + A_Minus(U.vars(i + 1, 0), grid.RNorms(i, 0)) * U.vars(i + 1, 0)) * grid.RArea(i, 0)			// Right Side
		+ (A_Plus(U.vars(i, 0), grid.BNorms(i, 0)) * U.vars(i, 0) + A_Minus(bottomBC.first, grid.BNorms(i, 0)) * bottomBC.first) * grid.BArea(i, 0)				// Bottom Side
		+ (A_Plus(U.vars(i, 0), grid.TNorms(i, 0)) * U.vars(i, 0) + A_Minus(U.vars(i, 1), grid.TNorms(i, 0)) * U.vars(i, 1)) * grid.TArea(i, 0)					// Top boundary 
		+ A_Minus(U.vars(i - 1, 0), grid.LNorms(i, 0)) * dU_old.vars(i - 1, 0) * grid.LArea(i, 0)								  								// Left dU addition
		+ A_Minus(U.vars(i + 1, 0), grid.RNorms(i, 0)) * dU_old.vars(i + 1, 0) * grid.RArea(i, 0)								 								// Right dU addition
		);


	alpha = A - B * g.strip(1);
	v.strip(0) = (F - B * v.strip(1)) / alpha;

	// Calculate dU
	dU.vars(0) = v.strip(0);

	for (int j = 1; j < Ny; ++j) {
		dU.vars(j) = v.strip(j) - g.strip(j) * dU.vars(j - 1);
	}

	return dU;
}


Matrix1D solveLeftLine(Matrix2D& U, Matrix U_inlet, Matrix2D& dU_old, const Grid& grid, 
	BoundaryCondition leftType, BoundaryCondition toptype, BoundaryCondition bottomtype, 
	const int i, const int Nx, const int Ny, const double dt) {

	Matrix A; // Relates to U(j+1) 
	Matrix B; // Relates to U(j) 
	Matrix C; // Relates to U(j-1) 
	Matrix F; // Right hand side 

	// Intermediate matrices
	Matrix alpha;
	Matrix2D v(4, 1, Nx);
	Matrix2D g(4, 4, Nx);

	Matrix1D dU(4, Ny);

	// Grab boundary values
	Duo bottomBC = Boundary2D(bottomtype, U.vars(i, 0), U_inlet, grid.BNorms(i, 0));  
	Duo topBC = Boundary2D(toptype, U.vars(i, Ny - 1), U_inlet, grid.BNorms(i, Ny - 1));


	// Top Boundary

	Duo leftBC = Boundary2D(leftType, U.vars(i, Ny - 1), U_inlet, grid.LNorms(i, Ny - 1));

	A = grid.Volume(i, Ny - 1) / dt * identity(4)
		+ (A_Plus(U.vars(i, Ny - 1), grid.TNorms(i, Ny - 1)) + topBC.second * A_Minus(topBC.first, grid.TNorms(i, Ny - 1))) * grid.TArea(i, Ny - 1)
		+ (A_Plus(U.vars(i, Ny - 1), grid.LNorms(i, Ny - 1)) + leftBC.second * A_Minus(leftBC.first, grid.LNorms(i, Ny - 1))) * grid.LArea(i, Ny - 1) 
		+ A_Plus(U.vars(i, Ny - 1), grid.RNorms(i, Ny - 1)) * grid.RArea(i, Ny - 1)
		+ A_Plus(U.vars(i, Ny - 1), grid.BNorms(i, Ny - 1)) * grid.BArea(i, Ny - 1);
		

	C = A_Minus(U.vars(i, Ny - 2), grid.BNorms(i, Ny - 1)) * grid.BArea(i, Ny - 1);

	F = -1 * ((A_Plus(U.vars(i, Ny - 1), grid.LNorms(i, Ny - 1)) * U.vars(i, Ny - 1) + A_Minus(leftBC.first, grid.LNorms(i, Ny - 1)) * leftBC.first) * grid.LArea(i, Ny - 1)				// Left Side		 
		+ (A_Plus(U.vars(i, Ny - 1), grid.RNorms(i, Ny - 1)) * U.vars(i, Ny - 1) + A_Minus(U.vars(i + 1, Ny - 1), grid.RNorms(i, Ny - 1)) * U.vars(i + 1, Ny - 1)) * grid.RArea(i, Ny - 1)	// Right Side
		+ (A_Plus(U.vars(i, Ny - 1), grid.BNorms(i, Ny - 1)) * U.vars(i, Ny - 1) + A_Minus(U.vars(i, Ny - 2), grid.BNorms(i, Ny - 1)) * U.vars(i, Ny - 2)) * grid.BArea(i, Ny - 1)			// Bottom Side
		+ (A_Plus(U.vars(i, Ny - 1), grid.TNorms(i, Ny - 1)) * U.vars(i, Ny - 1) + A_Minus(topBC.first, grid.TNorms(i, Ny - 1)) * topBC.first) * grid.TArea(i, Ny - 1)						// Top boundary
		+ A_Minus(U.vars(i + 1, Ny - 1), grid.RNorms(i, Ny - 1)) * dU_old.vars(i + 1, Ny - 1) * grid.RArea(i, Ny - 1)																		// Right dU addition
		);

	alpha = A;
	v.strip(Ny - 1) = F / alpha;
	g.strip(Ny - 1) = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		leftBC = Boundary2D(leftType, U.vars(i, j), U_inlet, grid.LNorms(i, j)); 

		B = A_Minus(U.vars(i, j + 1), grid.TNorms(i, j)) * grid.TArea(i, j);

		A = grid.Volume(i, j) / dt * identity(4)
			+ (A_Plus(U.vars(i, j), grid.LNorms(i, j)) + leftBC.second * A_Minus(leftBC.first, grid.LNorms(i, j))) * grid.LArea(i, j)  
			+ A_Plus(U.vars(i, j), grid.RNorms(i, j)) * grid.RArea(i, j)
			+ A_Plus(U.vars(i, j), grid.BNorms(i, j)) * grid.BArea(i, j)
			+ A_Plus(U.vars(i, j), grid.TNorms(i, j)) * grid.TArea(i, j);

		C = A_Minus(U.vars(i, j - 1), grid.BNorms(i, j)) * grid.BArea(i, j);

		F = -1 * ( (A_Plus(U.vars(i, j), grid.LNorms(i, j)) * U.vars(i, j) + A_Minus(leftBC.first, grid.LNorms(i, j)) * leftBC.first) * grid.LArea(i, j)	// Left Side	  	   
			+ (A_Plus(U.vars(i, j), grid.RNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i + 1, j), grid.RNorms(i, j)) * U.vars(i + 1, j)) * grid.RArea(i, j)	// Right Side
			+ (A_Plus(U.vars(i, j), grid.BNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i, j - 1), grid.BNorms(i, j)) * U.vars(i, j - 1)) * grid.BArea(i, j)	// Bottom Side
			+ (A_Plus(U.vars(i, j), grid.TNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i, j + 1), grid.TNorms(i, j)) * U.vars(i, j + 1)) * grid.TArea(i, j) // Top boundary 
			+ A_Minus(U.vars(i + 1, j), grid.RNorms(i, j)) * dU_old.vars(i + 1, j) * grid.RArea(i, j)														// Right dU addition
			);

		alpha = A - B * g.strip(j + 1);
		g.strip(j) = C / alpha;
		v.strip(j) = (F - B * v.strip(j + 1)) / alpha;
	}

	//  Bottom boundary
	leftBC = Boundary2D(leftType, U.vars(i, 0), U_inlet, grid.LNorms(i, 0)); 

	B = A_Minus(U.vars(i, 1), grid.TNorms(i, 0)) * grid.TArea(i, 0);

	A = grid.Volume(i, 0) / dt * identity(4)
		+ (A_Plus(U.vars(i, 0), grid.LNorms(i, 0)) + leftBC.second * A_Minus(leftBC.first, grid.LNorms(i, 0))) * grid.LArea(i, 0)  
		+ A_Plus(U.vars(i, 0), grid.RNorms(i, 0)) * grid.RArea(i, 0)
		+ A_Plus(U.vars(i, 0), grid.BNorms(i, 0)) * grid.BArea(i, 0)
		+ A_Plus(U.vars(i, 0), grid.TNorms(i, 0)) * grid.TArea(i, 0);

	F = -1 * ( (A_Plus(U.vars(i, 0), grid.LNorms(i, 0)) * U.vars(i, 0) + A_Minus(leftBC.first, grid.LNorms(i, 0)) * leftBC.first) * grid.LArea(i, 0)		// Left Side	   	   
		+ (A_Plus(U.vars(i, 0), grid.RNorms(i, 0)) * U.vars(i, 0) + A_Minus(U.vars(i + 1, 0), grid.RNorms(i, 0)) * U.vars(i + 1, 0)) * grid.RArea(i, 0)		// Right Side
		+ (A_Plus(U.vars(i, 0), grid.BNorms(i, 0)) * U.vars(i, 0) + A_Minus(bottomBC.first, grid.BNorms(i, 0)) * bottomBC.first) * grid.BArea(i, 0)			// Bottom Side
		+ (A_Plus(U.vars(i, 0), grid.TNorms(i, 0)) * U.vars(i, 0) + A_Minus(U.vars(i, 1), grid.TNorms(i, 0)) * U.vars(i, 1)) * grid.TArea(i, 0)				// Top boundary 
		+ A_Minus(U.vars(i + 1, 0), grid.RNorms(i, 0)) * dU_old.vars(i + 1, 0) * grid.RArea(i, 0)								 							// Right dU addition
		);


	alpha = A - B * g.strip(1);
	v.strip(0) = (F - B * v.strip(1)) / alpha;

	// Calculate dU
	dU.vars(0) = v.strip(0);

	for (int j = 1; j < Ny; ++j) {
		dU.vars(j) = v.strip(j) - g.strip(j) * dU.vars(j - 1);
	}

	return dU;
}

Matrix1D solveRightLine(Matrix2D& U, Matrix U_inlet, Matrix2D& dU_old, const Grid& grid, 
	BoundaryCondition righttype, BoundaryCondition toptype, BoundaryCondition bottomtype, 
	const int i, const int Nx, const int Ny, const double dt) {

	Matrix A; // Relates to U(j+1) 
	Matrix B; // Relates to U(j) 
	Matrix C; // Relates to U(j-1) 
	Matrix F; // Right hand side 

	// Intermediate matrices
	Matrix alpha;
	Matrix2D v(4, 1, Nx);
	Matrix2D g(4, 4, Nx);

	Matrix1D dU(4, Ny);

	// Grab boundary values
	Duo bottomBC = Boundary2D(bottomtype, U.vars(i, 0), U_inlet, grid.BNorms(i, 0));
	Duo topBC = Boundary2D(toptype, U.vars(i, Ny - 1), U_inlet, grid.TNorms(i, Ny - 1));


	// Top Boundary

	Duo rightBC = Boundary2D(righttype, U.vars(i, Ny - 1), U_inlet, grid.RNorms(i, Ny - 1));

	A = grid.Volume(i, Ny - 1) / dt * identity(4)  
		+ (A_Plus(U.vars(i, Ny - 1), grid.TNorms(i, Ny - 1)) + topBC.second * A_Minus(topBC.first, grid.TNorms(i, Ny - 1))) * grid.TArea(i, Ny - 1) 
		+ (A_Plus(U.vars(i, Ny - 1), grid.RNorms(i, Ny - 1)) + rightBC.second * A_Minus(rightBC.first, grid.RNorms(i, Ny - 1))) * grid.RArea(i, Ny - 1) 
		+ A_Plus(U.vars(i, Ny - 1), grid.BNorms(i, Ny - 1)) * grid.BArea(i, Ny - 1)
		+ A_Plus(U.vars(i, Ny - 1), grid.LNorms(i, Ny - 1)) * grid.LArea(i, Ny - 1);

	C = A_Minus(U.vars(i, Ny - 2), grid.BNorms(i, Ny - 1)) * grid.BArea(i, Ny - 1);

	F = -1 * ((A_Plus(U.vars(i, Ny - 1), grid.LNorms(i, Ny - 1)) * U.vars(i, Ny - 1) + A_Minus(U.vars(i - 1, Ny - 1), grid.LNorms(i, Ny - 1)) * U.vars(i - 1, Ny - 1)) * grid.LArea(i, Ny - 1)	// Left Side		
		+ (A_Plus(U.vars(i, Ny - 1), grid.RNorms(i, Ny - 1)) * U.vars(i, Ny - 1) + A_Minus(rightBC.first, grid.RNorms(i, Ny - 1)) * rightBC.first) * grid.RArea(i, Ny - 1)						// Right Side 
		+ (A_Plus(U.vars(i, Ny - 1), grid.BNorms(i, Ny - 1)) * U.vars(i, Ny - 1) + A_Minus(U.vars(i, Ny - 2), grid.BNorms(i, Ny - 1)) * U.vars(i, Ny - 2)) * grid.BArea(i, Ny - 1)				// Bottom Side
		+ (A_Plus(U.vars(i, Ny - 1), grid.TNorms(i, Ny - 1)) * U.vars(i, Ny - 1) + A_Minus(topBC.first, grid.TNorms(i, Ny - 1)) * topBC.first) * grid.TArea(i, Ny - 1)							// Top boundary
		+ A_Minus(U.vars(i - 1, Ny - 1), grid.LNorms(i, Ny - 1)) * dU_old.vars(i - 1, Ny - 1) * grid.LArea(i, Ny - 1)																			// Left dU addition
		);

	alpha = A;
	v.strip(Ny - 1) = F / alpha;
	g.strip(Ny - 1) = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		rightBC = Boundary2D(righttype, U.vars(i, j), U_inlet, grid.RNorms(i, j));

		B = A_Minus(U.vars(i, j + 1), grid.TNorms(i, j)) * grid.TArea(i, j); 

		A = grid.Volume(i, j) / dt * identity(4)
			+ A_Plus(U.vars(i, j), grid.LNorms(i, j)) * grid.LArea(i, j)
			+ (A_Plus(U.vars(i, j), grid.RNorms(i, j)) + rightBC.second * A_Minus(rightBC.first, grid.RNorms(i, j))) * grid.RArea(i, j)  
			+ A_Plus(U.vars(i, j), grid.BNorms(i, j)) * grid.BArea(i, j)
			+ A_Plus(U.vars(i, j), grid.TNorms(i, j)) * grid.TArea(i, j);

		C = A_Minus(U.vars(i, j - 1), grid.BNorms(i, j)) * grid.BArea(i, j);

		F = -1 * ((A_Plus(U.vars(i, j), grid.LNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i - 1, j), grid.LNorms(i, j)) * U.vars(i - 1, j)) * grid.LArea(i, j)		// Left Side		
			+ (A_Plus(U.vars(i, j), grid.RNorms(i, j)) * U.vars(i, j) + A_Minus(rightBC.first, grid.RNorms(i, j)) * rightBC.first) * grid.RArea(i, j)				// Right Side   
			+ (A_Plus(U.vars(i, j), grid.BNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i, j - 1), grid.BNorms(i, j)) * U.vars(i, j - 1)) * grid.BArea(i, j)			// Bottom Side
			+ (A_Plus(U.vars(i, j), grid.TNorms(i, j)) * U.vars(i, j) + A_Minus(U.vars(i, j + 1), grid.TNorms(i, j)) * U.vars(i, j + 1)) * grid.TArea(i, j)			// Top boundary 
			+ A_Minus(U.vars(i - 1, j), grid.LNorms(i, j)) * dU_old.vars(i - 1, j) * grid.LArea(i, j)																// Left dU addition
			);

		alpha = A - B * g.strip(j + 1);
		g.strip(j) = C / alpha;
		v.strip(j) = (F - B * v.strip(j + 1)) / alpha;
	}

	//  Bottom boundary
	rightBC = Boundary2D(righttype, U.vars(i, 0), U_inlet, grid.RNorms(i, 0)); 

	B = A_Minus(U.vars(i, 1), grid.TNorms(i, 0)) * grid.TArea(i, 0);

	A = grid.Volume(i, 0) / dt * identity(4)
		+ A_Plus(U.vars(i, 0), grid.LNorms(i, 0)) * grid.LArea(i, 0)
		+ (A_Plus(U.vars(i, 0), grid.RNorms(i, 0)) + rightBC.second * A_Minus(rightBC.first, grid.RNorms(i, 0))) * grid.RArea(i, 0) 
		+ A_Plus(U.vars(i, 0), grid.BNorms(i, 0)) * grid.BArea(i, 0)
		+ A_Plus(U.vars(i, 0), grid.TNorms(i, 0)) * grid.TArea(i, 0);

	F = -1 * ((A_Plus(U.vars(i, 0), grid.LNorms(i, 0)) * U.vars(i, 0) + A_Minus(U.vars(i - 1, 0), grid.LNorms(i, 0)) * U.vars(i - 1, 0)) * grid.LArea(i, 0)		// Left Side	 	
		+ (A_Plus(U.vars(i, 0), grid.RNorms(i, 0)) * U.vars(i, 0) + A_Minus(rightBC.first, grid.RNorms(i, 0)) * rightBC.first) * grid.RArea(i, 0)				// Right Side 
		+ (A_Plus(U.vars(i, 0), grid.BNorms(i, 0)) * U.vars(i, 0) + A_Minus(bottomBC.first, grid.BNorms(i, 0)) * bottomBC.first) * grid.BArea(i, 0)				// Bottom Side
		+ (A_Plus(U.vars(i, 0), grid.TNorms(i, 0)) * U.vars(i, 0) + A_Minus(U.vars(i, 1), grid.TNorms(i, 0)) * U.vars(i, 1)) * grid.TArea(i, 0)					// Top boundary 
		+ A_Minus(U.vars(i - 1, 0), grid.LNorms(i, 0)) * dU_old.vars(i - 1, 0) * grid.LArea(i, 0)								  								// Left dU addition
		);


	alpha = A - B * g.strip(1);
	v.strip(0) = (F - B * v.strip(1)) / alpha;

	// Calculate dU
	dU.vars(0) = v.strip(0);

	for (int j = 1; j < Ny; ++j) {
		dU.vars(j) = v.strip(j) - g.strip(j) * dU.vars(j - 1);
	}

	return dU;
}