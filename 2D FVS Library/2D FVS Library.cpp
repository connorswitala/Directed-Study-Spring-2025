// 2D FVS Library.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "framework.h"
#include "2DFVSLibrary.h"
#include <chrono>



// Function that converts primitives to conserved variables
Vector primtoCons(const Vector& V) {

	static Vector U(4, 0.0); 
	U[0] = V[0];
	U[1] = V[0] * V[1];
	U[2] = V[0] * V[2];
	U[3] = V[3] / (gamma - 1) + 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2]);

	return U;
}

// Function that converts conserved variables to primitives
Vector constoPrim(const Vector& U) {

	static Vector V(4, 0.0); 
	V[0] = U[0];
	V[1] = U[1] / U[0];
	V[2] = U[2] / U[0];
	V[3] = (U[3] - 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2])) * (gamma - 1);

	return V;
}

// Function that calculates ghost cells
Duo Boundary2D(BoundaryCondition type, const Vector& U, const Vector& U_inlet, const Point& normals) {   

	static Matrix E;
	static Vector Ug(4, 0.0), inside_Primitives(4, 0.0), ghost_Primitives(4, 0.0); 

	double u = 0.0, v = 0.0; 

	switch (type) {

	case BoundaryCondition::Inlet: 
		E = zeros(4, 4); 
		return make_pair(U_inlet, E);  

	case BoundaryCondition::Outlet: 
		E = ones(4, 4); 
		return make_pair(U, E); 

	case BoundaryCondition::Wall: 
		Ug = ones(4); 
		E = ones(4, 4); 		
		return make_pair(Ug, E); 

	case BoundaryCondition::Symmetry: 

		inside_Primitives = constoPrim(U); 
		u = inside_Primitives[1];
		v = inside_Primitives[2];	

		ghost_Primitives[0] = inside_Primitives[0]; 
		ghost_Primitives[1] = u - 2 * (u * normals.x + v * normals.y) * normals.x; 
		ghost_Primitives[2] = v - 2 * (u * normals.x + v * normals.y) * normals.y; 
		ghost_Primitives[3] = inside_Primitives[3]; 

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
	Tensor& U, 
	const Grid& grid,
	const int Nx,
	const int Ny) {

	Matrix density(Nx, vector<double>(Ny, 0.0)), u_vel(Nx, vector<double>(Ny, 0.0)), v_vel(Nx, vector<double>(Ny, 0.0)), pressure(Nx, vector<double>(Ny, 0.0)); 
	Vector Primitives(4); 

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) { 

			Primitives = constoPrim(U[i][j]); 
			density[i][j] = Primitives[0];
			u_vel[i][j] = Primitives[1];
			v_vel[i][j] = Primitives[2];
			pressure[i][j] = Primitives[3];
		}
	}

	ofstream file(filename);


	file << "density, u_velocity, v_velocity, pressure, x_points, y_points, z_points" << endl;

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			file << density[i][j] << ", " << u_vel[i][j] << ", " << v_vel[i][j] << ", " << pressure[i][j] << ", " << grid.Center(i, j).x << ", " << grid.Center(i, j).y << ", 0.0" << endl;
		}
	}

	file.close();
}

// Writes to CSV file (Line plot)
void writeCSV(const string& filename,
	const Tensor& U,
	const Grid& grid,
	const int j,
	const int Nx) {

	Vector density(Nx, 0.0), u_vel(Nx, 0.0), v_vel(Nx, 0.0), pressure(Nx, 0.0);   
	Vector Primitives(4);  

	for (int i = 0; i < Nx; ++i) {
			Primitives = constoPrim(U[i][j]); 
			density[i] = Primitives[0];
			u_vel[i] = Primitives[1];
			v_vel[i] = Primitives[2];
			pressure[i] = Primitives[3];
	}

	ofstream file(filename);

	file << "density, u_velocity, v_velocity, pressure, x_points, y_points" << endl;

	for (int i = 0; i < Nx; ++i) { 
		file << density[i] << ", " << u_vel[i] << ", " << v_vel[i] << ", " << pressure[i] << ", " << grid.Center(i, j).x << ", " << grid.Center(i, j).y << ", " << "0.0\n" << endl;  
	} 

	file.close();

}


void solveOneTimestep(Tensor& U, Tensor& dU_new, Vector U_inlet, Tensor& dU_old, const Grid& grid, 
	BoundaryConditions BoundaryTypes, const int& Nx, const int& Ny, double& dt, const double& CFL, 
	Tensor& i_Fluxes, Tensor& j_Fluxes, Tesseract& i_plus_Jacobians, Tesseract& i_minus_Jacobians, Tesseract& j_plus_Jacobians, Tesseract& j_minus_Jacobians) { 

	double inner_residual = 1.0;

	while (inner_residual >= 1e-8) {
		
		dt = calculate_dt(U, grid, Nx, Ny, CFL); 

		Calculate_Jacobians(U, U_inlet, BoundaryTypes, grid, Nx, Ny, i_Fluxes, j_Fluxes, i_plus_Jacobians, i_minus_Jacobians, j_plus_Jacobians, j_minus_Jacobians);

		solveLeftLine(U, dU_new, U_inlet, dU_old, grid, BoundaryTypes, 0, Nx, Ny, dt, i_Fluxes, j_Fluxes, i_plus_Jacobians, i_minus_Jacobians, j_plus_Jacobians, j_minus_Jacobians);

		for (int i = 1; i < Nx - 1; ++i) {
			solveMiddleLine(U, dU_new, U_inlet, dU_old, grid, BoundaryTypes, i, Nx, Ny, dt, i_Fluxes, j_Fluxes, i_plus_Jacobians, i_minus_Jacobians, j_plus_Jacobians, j_minus_Jacobians);
		}

		solveRightLine(U, dU_new, U_inlet, dU_old, grid, BoundaryTypes, Nx - 1, Nx, Ny, dt, i_Fluxes, j_Fluxes, i_plus_Jacobians, i_minus_Jacobians, j_plus_Jacobians, j_minus_Jacobians);

		inner_residual = calculateInnerResidual(U, dU_new, dU_old, grid, Nx, Ny, dt, i_Fluxes, j_Fluxes, i_plus_Jacobians, i_minus_Jacobians, j_plus_Jacobians, j_minus_Jacobians); 

		dU_old = dU_new;
	}

	U = U + dU_new; 

}

 //Function that calculates dt
double calculate_dt(Tensor& U, const Grid& grid, const int& Nx, const int& Ny, const double& CFL) {

	Vector V(4);
	double dx, dy, c, dt_old, dt = 1;

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			V = constoPrim(U[i][j]);
			dx = min(grid.jArea(i, j), grid.jArea(i, j + 1));
			dy = min(grid.iArea(i, j), grid.iArea(i + 1, j));
			c = sqrt(gamma * V[3] / V[0]);

			dt_old = CFL / (fabs(V[1] / dx) + fabs(V[2] / dy) + c * sqrt(1 / (dx * dx) + 1 / (dy * dy)));

			if (dt_old < dt) dt = dt_old;

		}
	}
	return dt; 
}

void Calculate_Jacobians(const Tensor& U, Vector& U_inlet, const BoundaryConditions& BoundaryTypes, const Grid& grid, int Nx, int Ny,
	Tensor& i_Fluxes, Tensor& j_Fluxes, Tesseract& i_plus_Jacobians, Tesseract& i_minus_Jacobians, Tesseract& j_plus_Jacobians, Tesseract& j_minus_Jacobians) {

	Vector V(4), Vi(4), Vii(4), Up; 
	Matrix M(4, Vector(4)), dvdu(4, Vector(4)), dudv(4, Vector(4)); 
	double g = 5.72;
	double weight, dp, a, rho, u, v, p, nx, ny, uprime, l1, l2, l3, l4;


	// Calculate Jacobians and Explicit fluxes for i-faces on left boundary 
	for (int j = 0; j < Ny; ++j) {

		Duo leftBC = Boundary2D(BoundaryTypes.left, U[0][j], U_inlet, grid.iNorms(0, j)); 

		Vi = constoPrim(leftBC.first); 
		Vii = constoPrim(U[0][j]); 

		nx = grid.iNorms(0, j).x;  
		ny = grid.iNorms(0, j).y;  

		dp = fabs(Vii[3] - Vi[3]) / min(Vii[3], Vi[3]); 
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
		V = weight * Vi + (1 - weight) * Vii; 

		rho = V[0]; 
		u = V[1]; 
		v = V[2];
		p = V[3];

		a = sqrt(gamma * p / rho);
		
		uprime = u * nx + v * ny;
		l1 = 0.5 * (uprime - a + fabs(uprime - a));  
		l2 = 0.5 * (uprime + fabs(uprime)); 
		l3 = 0.5 * (uprime + fabs(uprime));
		l4 = 0.5 * (uprime + a + fabs(uprime + a)); 
		
		dvdu = {
		{1, 0, 0, 0},
		{-u / rho, 1 / rho, 0, 0},
		{-v / rho, 0, 1 / rho, 0},
		{0.5 * (gamma - 1) * (u * u + v * v), -u * (gamma - 1), -v * (gamma - 1), gamma - 1}
		};

		dudv = {
			{1, 0, 0, 0},
			{u, rho, 0, 0},
			{v, 0, rho, 0},
			{0.5 * (u * u + v * v), rho * u, rho * v, 1 / (gamma - 1)}
		};

		M = { 
		{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)}, 
		{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)}, 
		{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)}, 
		{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)} 
		}; 

		i_plus_Jacobians[0][j] = dudv * M * dvdu;

		l1 = 0.5 * (uprime - a - fabs(uprime - a));
		l2 = 0.5 * (uprime - fabs(uprime));
		l3 = 0.5 * (uprime - fabs(uprime));
		l4 = 0.5 * (uprime + a - fabs(uprime + a)); 

		M = {
		{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)},
		{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
		};

		i_minus_Jacobians[0][j] = dudv * M * dvdu;
		i_Fluxes[0][j] = i_plus_Jacobians[0][j] * leftBC.first + i_minus_Jacobians[0][j] * U[0][j]; 
	}

	// Calculate Jacobians and Explicit fluxes for i-faces on right boundary
	for (int j = 0; j < Ny; ++j) {

		Duo rightBC = Boundary2D(BoundaryTypes.right, U[Nx - 1][j], U_inlet, grid.iNorms(Nx, j));   

		Vi = constoPrim(U[Nx - 1][j]); 
		Vii = constoPrim(rightBC.first);

		nx = grid.iNorms(Nx, j).x; 
		ny = grid.iNorms(Nx, j).y;  

		dp = fabs(Vii[3] - Vi[3]) / min(Vii[3], Vi[3]);
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
		V = weight * Vi + (1 - weight) * Vii;

		rho = V[0];
		u = V[1];
		v = V[2];
		p = V[3];

		a = sqrt(gamma * p / rho);

		uprime = u * nx + v * ny;
		l1 = 0.5 * (uprime - a + fabs(uprime - a));
		l2 = 0.5 * (uprime + fabs(uprime));
		l3 = 0.5 * (uprime + fabs(uprime));
		l4 = 0.5 * (uprime + a + fabs(uprime + a));

		dvdu = {
		{1, 0, 0, 0},
		{-u / rho, 1 / rho, 0, 0},
		{-v / rho, 0, 1 / rho, 0},
		{0.5 * (gamma - 1) * (u * u + v * v), -u * (gamma - 1), -v * (gamma - 1), gamma - 1}
		};

		dudv = {
			{1, 0, 0, 0},
			{u, rho, 0, 0},
			{v, 0, rho, 0},
			{0.5 * (u * u + v * v), rho * u, rho * v, 1 / (gamma - 1)}
		};

		M = {
		{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)},
		{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
		};

		i_plus_Jacobians[Nx][j] = dudv * M * dvdu; 

		l1 = 0.5 * (uprime - a - fabs(uprime - a));
		l2 = 0.5 * (uprime - fabs(uprime));
		l3 = 0.5 * (uprime - fabs(uprime));
		l4 = 0.5 * (uprime + a - fabs(uprime + a));

		M = {
		{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)},
		{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
		};

		i_minus_Jacobians[Nx][j] = dudv * M * dvdu; 
		i_Fluxes[Nx][j] = i_plus_Jacobians[Nx][j] * U[Nx - 1][j] + i_minus_Jacobians[Nx][j] * rightBC.first; 
	}

	// Calculate Jacobians and Explicit fluxes for j-faces on bottom boundary
	for (int i = 0; i < Nx; ++i) { 

		Duo bottomBC = Boundary2D(BoundaryTypes.bottom, U[i][0], U_inlet, grid.jNorms(i, 0));  

		Vi = constoPrim(bottomBC.first); 
		Vii = constoPrim(U[i][0]); 

		nx = grid.jNorms(i, 0).x;  
		ny = grid.jNorms(i, 0).y;  

		dp = fabs(Vii[3] - Vi[3]) / min(Vii[3], Vi[3]);
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
		V = weight * Vi + (1 - weight) * Vii;

		rho = V[0];
		u = V[1];
		v = V[2];
		p = V[3];

		a = sqrt(gamma * p / rho);

		uprime = u * nx + v * ny;
		l1 = 0.5 * (uprime - a + fabs(uprime - a));
		l2 = 0.5 * (uprime + fabs(uprime));
		l3 = 0.5 * (uprime + fabs(uprime));
		l4 = 0.5 * (uprime + a + fabs(uprime + a));

		dvdu = {
		{1, 0, 0, 0},
		{-u / rho, 1 / rho, 0, 0},
		{-v / rho, 0, 1 / rho, 0},
		{0.5 * (gamma - 1) * (u * u + v * v), -u * (gamma - 1), -v * (gamma - 1), gamma - 1}
		};

		dudv = {
			{1, 0, 0, 0},
			{u, rho, 0, 0},
			{v, 0, rho, 0},
			{0.5 * (u * u + v * v), rho * u, rho * v, 1 / (gamma - 1)}
		};

		M = {
		{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)},
		{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
		};

		j_plus_Jacobians[i][0] = dudv * M * dvdu; 

		l1 = 0.5 * (uprime - a - fabs(uprime - a));
		l2 = 0.5 * (uprime - fabs(uprime));
		l3 = 0.5 * (uprime - fabs(uprime));
		l4 = 0.5 * (uprime + a - fabs(uprime + a));

		M = {
		{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)},
		{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
		};

		j_minus_Jacobians[i][0] = dudv * M * dvdu; 
		j_Fluxes[i][0] = j_plus_Jacobians[i][0] * bottomBC.first + j_minus_Jacobians[i][0] * U[i][0]; 

	}

	// Calculate Jacobians and Explicit fluxes for j-faces on top boundary
	for (int i = 0; i < Nx; ++i) {

		Duo topBC = Boundary2D(BoundaryTypes.top, U[i][Ny - 1], U_inlet, grid.jNorms(i, Ny)); 

		Vi = constoPrim(topBC.first);
		Vii = constoPrim(U[i][Ny - 1]); 

		nx = grid.jNorms(i, Ny).x;   
		ny = grid.jNorms(i, Ny).y;   


		dp = fabs(Vii[3] - Vi[3]) / min(Vii[3], Vi[3]);
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
		V = weight * Vi + (1 - weight) * Vii;

		rho = V[0];
		u = V[1];
		v = V[2];
		p = V[3];

		a = sqrt(gamma * p / rho);

		uprime = u * nx + v * ny;
		l1 = 0.5 * (uprime - a + fabs(uprime - a));
		l2 = 0.5 * (uprime + fabs(uprime));
		l3 = 0.5 * (uprime + fabs(uprime));
		l4 = 0.5 * (uprime + a + fabs(uprime + a));

		dvdu = {
		{1, 0, 0, 0},
		{-u / rho, 1 / rho, 0, 0},
		{-v / rho, 0, 1 / rho, 0},
		{0.5 * (gamma - 1) * (u * u + v * v), -u * (gamma - 1), -v * (gamma - 1), gamma - 1}
		};

		dudv = {
			{1, 0, 0, 0},
			{u, rho, 0, 0},
			{v, 0, rho, 0},
			{0.5 * (u * u + v * v), rho * u, rho * v, 1 / (gamma - 1)}
		};

		M = {
		{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)},
		{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
		};

		j_plus_Jacobians[i][Ny] = dudv * M * dvdu; 

		l1 = 0.5 * (uprime - a - fabs(uprime - a));
		l2 = 0.5 * (uprime - fabs(uprime));
		l3 = 0.5 * (uprime - fabs(uprime));
		l4 = 0.5 * (uprime + a - fabs(uprime + a));

		M = {
		{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)},
		{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)},
		{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
		};

		j_minus_Jacobians[i][Ny] = dudv * M * dvdu; 
		j_Fluxes[i][Ny] = j_plus_Jacobians[i][Ny] * U[i][Ny - 1] + j_minus_Jacobians[i][Ny] * topBC.first; 

	}

	// Inner i-faces  
	for (int i = 1; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {

			Vi = constoPrim(U[i - 1][j]);
			Vii = constoPrim(U[i][j]);

			nx = grid.iNorms(i, j).x;
			ny = grid.iNorms(i, j).y;

			dp = fabs(Vii[3] - Vi[3]) / min(Vii[3], Vi[3]);
			weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
			V = weight * Vi + (1 - weight) * Vii;

			rho = V[0];
			u = V[1];
			v = V[2];
			p = V[3];

			a = sqrt(gamma * p / rho);

			uprime = u * nx + v * ny;
			l1 = 0.5 * (uprime - a + fabs(uprime - a));
			l2 = 0.5 * (uprime + fabs(uprime));
			l3 = 0.5 * (uprime + fabs(uprime));
			l4 = 0.5 * (uprime + a + fabs(uprime + a));

			dvdu = {
			{1, 0, 0, 0},
			{-u / rho, 1 / rho, 0, 0},
			{-v / rho, 0, 1 / rho, 0},
			{0.5 * (gamma - 1) * (u * u + v * v), -u * (gamma - 1), -v * (gamma - 1), gamma - 1}
			};

			dudv = {
				{1, 0, 0, 0},
				{u, rho, 0, 0},
				{v, 0, rho, 0},
				{0.5 * (u * u + v * v), rho * u, rho * v, 1 / (gamma - 1)}
			};

			M = {
			{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)},
			{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)},
			{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)},
			{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
			};

			i_plus_Jacobians[i][j] = dudv * M * dvdu;

			l1 = 0.5 * (uprime - a - fabs(uprime - a));
			l2 = 0.5 * (uprime - fabs(uprime));
			l3 = 0.5 * (uprime - fabs(uprime));
			l4 = 0.5 * (uprime + a - fabs(uprime + a));

			M = {
			{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)},
			{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)},
			{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)},
			{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
			};

			i_minus_Jacobians[i][j] = dudv * M * dvdu;
			i_Fluxes[i][j] = i_plus_Jacobians[i][j] * U[i - 1][j] + i_minus_Jacobians[i][j] * U[i][j];

		}
	}

	// Inner j-faces 
	for (int i = 0; i < Nx; ++i) {
		for (int j = 1; j < Ny; ++j) {
			
			Vi = constoPrim(U[i][j - 1]); 
			Vii = constoPrim(U[i][j]); 

			nx = grid.jNorms(i, j).x;  
			ny = grid.jNorms(i, j).y;  

			dp = fabs(Vii[3] - Vi[3]) / min(Vii[3], Vi[3]);
			weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
			V = weight * Vi + (1 - weight) * Vii;

			rho = V[0];
			u = V[1];
			v = V[2];
			p = V[3];

			a = sqrt(gamma * p / rho);
			uprime = u * nx + v * ny;

			l1 = 0.5 * (uprime - a + fabs(uprime - a));
			l2 = 0.5 * (uprime + fabs(uprime));
			l3 = 0.5 * (uprime + fabs(uprime));
			l4 = 0.5 * (uprime + a + fabs(uprime + a));

			dvdu = {
			{1, 0, 0, 0},
			{-u / rho, 1 / rho, 0, 0},
			{-v / rho, 0, 1 / rho, 0},
			{0.5 * (gamma - 1) * (u * u + v * v), -u * (gamma - 1), -v * (gamma - 1), gamma - 1}
			};

			dudv = {
				{1, 0, 0, 0},
				{u, rho, 0, 0},
				{v, 0, rho, 0},
				{0.5 * (u * u + v * v), rho * u, rho * v, 1 / (gamma - 1)}
			};

			M = {
			{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)},
			{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)},
			{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)},
			{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
			};

			j_plus_Jacobians[i][j] = dudv * M * dvdu; 

			l1 = 0.5 * (uprime - a - fabs(uprime - a));
			l2 = 0.5 * (uprime - fabs(uprime));
			l3 = 0.5 * (uprime - fabs(uprime));
			l4 = 0.5 * (uprime + a - fabs(uprime + a));

			M = {
			{l2, rho * nx / (2 * a) * (-l1 + l4), rho * ny / (2 * a) * (-l1 + l4), 1 / (2 * a * a) * (l1 - 2 * l2 + l4)},
			{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * a * rho) * (-l1 + l4)},
			{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * a * rho) * (-l1 + l4)},
			{0, 0.5 * a * rho * nx * (-l1 + l4), 0.5 * a * rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
			};

			j_minus_Jacobians[i][j] = dudv * M * dvdu;  
			j_Fluxes[i][j] = j_plus_Jacobians[i][j] * U[i][j - 1] + j_minus_Jacobians[i][j] * U[i][j]; 
		
		}
	}

}

void solveLeftLine(Tensor& U, Tensor& dU_new, Vector U_inlet, Tensor& dU_old, const Grid& grid,
	BoundaryConditions BoundaryType, const int i, const int Nx, const int Ny, const double dt,
	Tensor& i_Fluxes, Tensor& j_Fluxes, Tesseract& i_plus_Jacobians, Tesseract& i_minus_Jacobians, Tesseract& j_plus_Jacobians, Tesseract& j_minus_Jacobians) {

	static Matrix A; // Relates to U(j+1) 
	static Matrix B; // Relates to U(j) 
	static Matrix C; // Relates to U(j-1) 
	static Vector F; // Right hand side  

	// Intermediate matrices
	static Matrix alpha;
	static Matrix v(Ny, Vector(4));
	static Tensor g(Ny, Matrix(4, Vector(4)));

	// Grab boundary values
	Duo bottomBC = Boundary2D(BoundaryType.bottom, U[i][0], U_inlet, grid.jNorms(i, 0));
	Duo topBC = Boundary2D(BoundaryType.top, U[i][Ny - 1], U_inlet, grid.jNorms(i, Ny));

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity(4)
		- i_minus_Jacobians[i][Ny - 1] * grid.iArea(i, Ny - 1)
		+ i_plus_Jacobians[i + 1][Ny - 1] * grid.iArea(i + 1, Ny - 1)
		- j_minus_Jacobians[i][Ny - 1] * grid.jArea(i, Ny - 1)
		+ (j_plus_Jacobians[i][Ny] + topBC.second * j_minus_Jacobians[i][Ny]) * grid.jArea(i, Ny);

	C = j_plus_Jacobians[i][Ny - 1] * (-grid.jArea(i, Ny - 1));

	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1))
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny))
		+ i_minus_Jacobians[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1)) * dU_old[i + 1][Ny - 1];

	alpha = A;
	v[Ny - 1] = F / alpha;
	g[Ny - 1] = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		B = j_minus_Jacobians[i][j + 1] * (grid.jArea(i, j + 1));

		A = grid.Volume(i, j) / dt * identity(4)
			- i_minus_Jacobians[i][j] * grid.iArea(i, j)
			+ i_plus_Jacobians[i + 1][j] * grid.iArea(i + 1, j)
			- j_minus_Jacobians[i][j] * grid.jArea(i, j)
			+ j_plus_Jacobians[i][j + 1] * grid.jArea(i, j + 1);

		C = j_plus_Jacobians[i][j] * (-grid.jArea(i, j));

		F = i_Fluxes[i][j] * (grid.iArea(i, j))
			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j))
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1))
			+ i_minus_Jacobians[i + 1][j] * (-grid.iArea(i + 1, j)) * dU_old[i + 1][j];



		alpha = A - B * g[j + 1];
		g[j] = C / alpha;
		v[j] = (F - B * v[j + 1]) / alpha;
	}

	//  Bottom boundary
	B = j_minus_Jacobians[i][1] * (-grid.jArea(i, 1));

	A = grid.Volume(i, 0) / dt * identity(4)
		- i_minus_Jacobians[i][0] * grid.iArea(i, 0)
		+ i_plus_Jacobians[i + 1][0] * grid.iArea(i + 1, 0)
		- (bottomBC.second * j_plus_Jacobians[i][0] + j_minus_Jacobians[i][0]) * grid.jArea(i, 0)
		+ j_plus_Jacobians[i][1] * grid.jArea(i, 1);

	F = i_Fluxes[i][0] * (grid.iArea(i, 0))
		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1))
		+ i_minus_Jacobians[i + 1][0] * (-grid.iArea(i + 1, 0)) * dU_old[i + 1][0];


	alpha = A - B * g[1];
	v[0] = (F - B * v[1]) / alpha;

	// Calculate dU
	dU_new[i][0] = v[0];

	for (int j = 1; j < Ny; ++j) {
		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];
	}

}

void solveMiddleLine(Tensor& U, Tensor& dU_new, Vector U_inlet, Tensor& dU_old, const Grid& grid,   
	BoundaryConditions BoundaryType, const int i, const int Nx, const int Ny, const double dt,
	Tensor& i_Fluxes, Tensor& j_Fluxes, Tesseract& i_plus_Jacobians, Tesseract& i_minus_Jacobians, Tesseract& j_plus_Jacobians, Tesseract& j_minus_Jacobians) {

	//auto line1 = chrono::high_resolution_clock::now(); // Start the timer  

	static Matrix A; // Relates to U(j+1) 
	static Matrix B; // Relates to U(j) 
	static Matrix C; // Relates to U(j-1) 
	static Vector F; // Right hand side  

	// Intermediate matrices
	static Matrix alpha; 
	static Matrix v(Ny, Vector(4));  
	static Tensor g(Ny, Matrix(4, Vector(4)));  

	// Grab boundary values
	Duo bottomBC = Boundary2D(BoundaryType.bottom, U[i][0], U_inlet, grid.jNorms(i, 0)); 
	Duo topBC = Boundary2D(BoundaryType.top, U[i][Ny - 1], U_inlet, grid.jNorms(i, Ny));   

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity(4) 
		- i_minus_Jacobians[i][Ny - 1] * grid.iArea(i, Ny - 1) 
		+ i_plus_Jacobians[i + 1][Ny - 1] * grid.iArea(i + 1, Ny - 1)
		- j_minus_Jacobians[i][Ny - 1] * grid.jArea(i, Ny - 1) 
		+ (j_plus_Jacobians[i][Ny] + topBC.second * j_minus_Jacobians[i][Ny]) * grid.jArea(i, Ny); 

	C = j_plus_Jacobians[i][Ny - 1] * (-grid.jArea(i, Ny - 1)); 

	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1)) 
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny))
		+ i_plus_Jacobians[i][Ny - 1] * (grid.iArea(i, Ny - 1)) * dU_old[i - 1][Ny - 1]
		+ i_minus_Jacobians[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1)) * dU_old[i + 1][Ny - 1]; 
	
	alpha = A;
	v[Ny - 1] = F / alpha;  
	g[Ny - 1] = C / alpha; 

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		B = j_minus_Jacobians[i][j + 1] * (grid.jArea(i, j + 1)); 

		A = grid.Volume(i, j) / dt * identity(4) 
			- i_minus_Jacobians[i][j] * grid.iArea(i, j) 
			+ i_plus_Jacobians[i + 1][j] * grid.iArea(i + 1, j)  
			- j_minus_Jacobians[i][j] * grid.jArea(i, j)  
			+ j_plus_Jacobians[i][j + 1] * grid.jArea(i, j + 1);   

		C = j_plus_Jacobians[i][j] * (-grid.jArea(i, j));  

		F = i_Fluxes[i][j] * (grid.iArea(i, j))  
			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j)) 
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1)) 
			+ i_plus_Jacobians[i][j] * (grid.iArea(i, j)) * dU_old[i - 1][j] 
			+ i_minus_Jacobians[i + 1][j] * (-grid.iArea(i + 1, j)) * dU_old[i + 1][j]; 



		alpha = A - B * g[j + 1]; 
		g[j] = C / alpha; 
		v[j] = (F - B * v[j + 1]) / alpha; 
	}

	//  Bottom boundary
	B = j_minus_Jacobians[i][1] * (-grid.jArea(i, 1));

	A = grid.Volume(i, 0) / dt * identity(4)
		- i_minus_Jacobians[i][0] * grid.iArea(i, 0)
		+ i_plus_Jacobians[i + 1][0] * grid.iArea(i + 1, 0)
		- (bottomBC.second * j_plus_Jacobians[i][0] + j_minus_Jacobians[i][0]) * grid.jArea(i, 0) 
		+ j_plus_Jacobians[i][1] * grid.jArea(i, 1); 

	F = i_Fluxes[i][0] * (grid.iArea(i, 0))
		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1))
		+ i_plus_Jacobians[i][0] * (grid.iArea(i, 0)) * dU_old[i - 1][0]
		+ i_minus_Jacobians[i + 1][0] * (-grid.iArea(i + 1, 0)) * dU_old[i + 1][0];


	alpha = A - B * g[1];
	v[0] = (F - B * v[1]) / alpha; 

	// Calculate dU
	dU_new[i][0] = v[0]; 

	for (int j = 1; j < Ny; ++j) {
		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];  
	}

	//auto line2 = chrono::high_resolution_clock::now(); // Stop the timer 
	//chrono::duration<double> duration2 = line2 - line1; // Calculate the duration 
	//cout <<"Line " << i << " took " << duration2.count() << "seconds." << endl; 
}

void solveRightLine(Tensor& U, Tensor& dU_new, Vector U_inlet, Tensor& dU_old, const Grid& grid,
	BoundaryConditions BoundaryType, const int i, const int Nx, const int Ny, const double dt,
	Tensor& i_Fluxes, Tensor& j_Fluxes, Tesseract& i_plus_Jacobians, Tesseract& i_minus_Jacobians, Tesseract& j_plus_Jacobians, Tesseract& j_minus_Jacobians) {

	static Matrix A; // Relates to U(j+1) 
	static Matrix B; // Relates to U(j) 
	static Matrix C; // Relates to U(j-1) 
	static Vector F; // Right hand side  

	// Intermediate matrices
	static Matrix alpha;
	static Matrix v(Ny, Vector(4));
	static Tensor g(Ny, Matrix(4, Vector(4)));

	// Grab boundary values
	Duo bottomBC = Boundary2D(BoundaryType.bottom, U[i][0], U_inlet, grid.jNorms(i, 0));
	Duo topBC = Boundary2D(BoundaryType.top, U[i][Ny - 1], U_inlet, grid.jNorms(i, Ny));

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity(4)
		- i_minus_Jacobians[i][Ny - 1] * grid.iArea(i, Ny - 1)
		+ i_plus_Jacobians[i + 1][Ny - 1] * grid.iArea(i + 1, Ny - 1)
		- j_minus_Jacobians[i][Ny - 1] * grid.jArea(i, Ny - 1)
		+ (j_plus_Jacobians[i][Ny] + topBC.second * j_minus_Jacobians[i][Ny]) * grid.jArea(i, Ny);

	C = j_plus_Jacobians[i][Ny - 1] * (-grid.jArea(i, Ny - 1));

	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1))
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny))
		+ i_plus_Jacobians[i][Ny - 1] * (grid.iArea(i, Ny - 1)) * dU_old[i - 1][Ny - 1];

	alpha = A;
	v[Ny - 1] = F / alpha;
	g[Ny - 1] = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		B = j_minus_Jacobians[i][j + 1] * (grid.jArea(i, j + 1));

		A = grid.Volume(i, j) / dt * identity(4)
			- i_minus_Jacobians[i][j] * grid.iArea(i, j)
			+ i_plus_Jacobians[i + 1][j] * grid.iArea(i + 1, j)
			- j_minus_Jacobians[i][j] * grid.jArea(i, j)
			+ j_plus_Jacobians[i][j + 1] * grid.jArea(i, j + 1);

		C = j_plus_Jacobians[i][j] * (-grid.jArea(i, j));

		F = i_Fluxes[i][j] * (grid.iArea(i, j))
			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j))
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1))
			+ i_plus_Jacobians[i][j] * (grid.iArea(i, j)) * dU_old[i - 1][j];



		alpha = A - B * g[j + 1];
		g[j] = C / alpha;
		v[j] = (F - B * v[j + 1]) / alpha;
	}

	//  Bottom boundary
	B = j_minus_Jacobians[i][1] * (-grid.jArea(i, 1));

	A = grid.Volume(i, 0) / dt * identity(4)
		- i_minus_Jacobians[i][0] * grid.iArea(i, 0)
		+ i_plus_Jacobians[i + 1][0] * grid.iArea(i + 1, 0)
		- (bottomBC.second * j_plus_Jacobians[i][0] + j_minus_Jacobians[i][0]) * grid.jArea(i, 0)
		+ j_plus_Jacobians[i][1] * grid.jArea(i, 1);

	F = i_Fluxes[i][0] * (grid.iArea(i, 0))
		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1))
		+ i_plus_Jacobians[i][0] * (grid.iArea(i, 0)) * dU_old[i - 1][0];


	alpha = A - B * g[1];
	v[0] = (F - B * v[1]) / alpha;

	// Calculate dU
	dU_new[i][0] = v[0];

	for (int j = 1; j < Ny; ++j) {
		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];
	}

}

double calculateInnerResidual(Tensor& U, Tensor& dU_new, Tensor& dU_old, const Grid& grid, const int Nx, const int Ny, const double dt,
	Tensor& i_Fluxes, Tensor& j_Fluxes, Tesseract& i_plus_Jacobians, Tesseract& i_minus_Jacobians, Tesseract& j_plus_Jacobians, Tesseract& j_minus_Jacobians) {
	//
	double inner_residual = 0.0;
	Vector res(4, 0.0);

	Matrix A; // Relates to U(j+1) 
	Matrix B; // Relates to U(j) 
	Matrix C; // Relates to U(j-1) 
	Vector F; // Right hand side  


	// Middle Boundary
	for (int i = 1; i < Nx - 1; ++i) {
		for (int j = 1; j < Ny - 1; ++j) {

			B = j_minus_Jacobians[i][j + 1] * (grid.jArea(i, j + 1));

			A = grid.Volume(i, j) / dt * identity(4)
				- i_minus_Jacobians[i][j] * grid.iArea(i, j)
				+ i_plus_Jacobians[i + 1][j] * grid.iArea(i + 1, j)
				- j_minus_Jacobians[i][j] * grid.jArea(i, j)
				+ j_plus_Jacobians[i][j + 1] * grid.jArea(i, j + 1);

			C = j_plus_Jacobians[i][j] * (-grid.jArea(i, j));

			F = i_Fluxes[i][j] * (grid.iArea(i, j))
				+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
				+ j_Fluxes[i][j] * (grid.jArea(i, j))
				+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1))
				+ i_plus_Jacobians[i][j] * (grid.iArea(i, j)) * dU_old[i - 1][j]
				+ i_minus_Jacobians[i + 1][j] * (-grid.iArea(i + 1, j)) * dU_old[i + 1][j];

			res = B * dU_new[i][j + 1] + A * dU_new[i][j] + C * dU_new[i][j - 1] - F;
			inner_residual = inner_residual + res[0] * res[0]; 

		}
	}
	return sqrt(inner_residual); 
}

double calculateResidual(const Tensor& U, const Grid& grid, int Nx, int Ny, const Tensor& i_Fluxes, const Tensor& j_Fluxes) {

	double res = 0.0;
	Vector intres(4, 0.0);

	for (int i = 1; i < Nx - 1; ++i) {
		for (int j = 1; j < Ny - 1; ++j) {
			intres = (i_Fluxes[i][j] * (-grid.iArea(i, j))			// Left Side		   
				+ i_Fluxes[i + 1][j] * (grid.iArea(i + 1, j))		// Right Side 
				+ j_Fluxes[i][j] * (-grid.jArea(i, j))			// Bottom Side  
				+ j_Fluxes[i][j + 1] * (grid.jArea(i, j + 1)))/ grid.Volume(i, j); 

			res = res + intres[0] * intres[0];
		
		}
	}

	return sqrt(res); 
}