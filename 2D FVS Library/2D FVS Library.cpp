// 2D FVS Library.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "framework.h"
#include "2DFVSLibrary.h"



// Function that converts primitives to conserved variables
Vector primtoCons(const Vector& V) {

	Vector U = zerosV(); 
	U[0] = V[0];
	U[1] = V[0] * V[1];
	U[2] = V[0] * V[2];
	U[3] = V[3] / (gamma - 1) + 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2]);

	return U;
}

// Function that converts conserved variables to primitives
Vector constoPrim(const Vector& U) {

	static Vector V = zerosV(); 
	V[0] = U[0];
	V[1] = U[1] / U[0];
	V[2] = U[2] / U[0];
	V[3] = (U[3] - 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2])) * (gamma - 1);

	return V;
}



// Function that calculates ghost cells
Duo Boundary2D(BoundaryCondition type, const Vector& U, const Vector& U_inlet, const Point& normals) {    

	Matrix E = zerosM();  
	Vector Ug = zerosV(), inside_Primitives = zerosV(), ghost_Primitives = zerosV(); 

	double u = 0.0, v = 0.0; 

	switch (type) {

	case BoundaryCondition::Inlet: 
		return make_pair(U_inlet, E);  

	case BoundaryCondition::Outlet: 
		E = onesM();  
		return make_pair(U, E); 

	case BoundaryCondition::Wall: 
		Ug = onesV(); 
		E = onesM(); 		 
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

		E = { {{1, 0, 0, 0},
			  {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
			  { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
			  { 0, 0, 0, 1 }}
		};

		return make_pair(Ug, E);

	default:
		throw invalid_argument("Unknown boundary condition type.");
	}

}


// Writes to VTK file (Contout plot)
void write2DCSV(const string& filename,  
	CellTensor& U, 
	const Grid& grid) {

	double a; 

	array<array<double, Ny>, Nx> density, u_velocity, v_velocity, pressure, Mach, Temperature; 
	Vector Primitives;

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) { 

			Primitives = constoPrim(U[i][j]); 
			density[i][j] = Primitives[0];
			pressure[i][j] = Primitives[3]; 
			a = sqrt(gamma * pressure[i][j] / density[i][j]);   
			u_velocity[i][j] = Primitives[1]; 
			v_velocity[i][j] = Primitives[2];
			Mach[i][j] = sqrt(u_velocity[i][j] * u_velocity[i][j] + v_velocity[i][j] * v_velocity[i][j]) / a; 
			Temperature[i][j] = pressure[i][j] / (density[i][j] * R);  
		}
	}

	ofstream file(filename);


	file << "density, u velocity, v velocity, pressure, Mach, Temperature, x_points, y_points, z_points" << endl; 

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			file << density[i][j] << ", " << u_velocity[i][j] << ", " << v_velocity[i][j] << ", " << pressure[i][j] << ", " << Mach[i][j] << ", " << Temperature[i][j] << ", " << grid.Center(i, j).x << ", " << grid.Center(i, j).y << ", 0.0" << endl;
		}
	}

	file.close();
	cout << "2D File saved successfully as \"" << filename << "\"" << endl; 
}

// Writes to CSV file (Line plot)
void write1DCSV(const string& filename,
	const CellTensor& U,
	const Grid& grid,
	const int j) {

	array<double, Nx> density, u_vel, v_vel, pressure; 
	Vector Primitives;  

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
	cout << "1D File written successfully." << endl;
}


void solveOneTimestep(CellTensor& U, CellTensor& dU_new, Vector U_inlet, CellTensor& dU_old, const Grid& grid,
	BoundaryConditions BoundaryTypes, double& dt, const double& CFL,
	iFaceTensor& i_Fluxes, jFaceTensor& j_Fluxes, iFaceTesseract& i_plus_Jacobians, iFaceTesseract& i_minus_Jacobians, jFaceTesseract& j_plus_Jacobians, jFaceTesseract& j_minus_Jacobians) {

	double inner_residual = 1.0;

	while (inner_residual >= 1e-8) {		

		Calculate_Jacobians(U, U_inlet, BoundaryTypes, grid, i_Fluxes, j_Fluxes, i_plus_Jacobians, i_minus_Jacobians, j_plus_Jacobians, j_minus_Jacobians);

		solveLeftLine(U, dU_new, U_inlet, dU_old, grid, BoundaryTypes, 0, dt, i_Fluxes, j_Fluxes, i_plus_Jacobians, i_minus_Jacobians, j_plus_Jacobians, j_minus_Jacobians);


		for (int i = 1; i < Nx - 1; ++i) {
			solveMiddleLine(U, dU_new, U_inlet, dU_old, grid, BoundaryTypes, i, dt, i_Fluxes, j_Fluxes, i_plus_Jacobians, i_minus_Jacobians, j_plus_Jacobians, j_minus_Jacobians);
		}

		solveRightLine(U, dU_new, U_inlet, dU_old, grid, BoundaryTypes, Nx - 1, dt, i_Fluxes, j_Fluxes, i_plus_Jacobians, i_minus_Jacobians, j_plus_Jacobians, j_minus_Jacobians);

	
		inner_residual = calculateInnerResidual(U, dU_new, dU_old, grid, dt, i_Fluxes, j_Fluxes, i_plus_Jacobians, i_minus_Jacobians, j_plus_Jacobians, j_minus_Jacobians); 
		dU_old = dU_new;
	}

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < 4; ++k) {
				U[i][j][k] += dU_new[i][j][k];   
			}
		}
	}

}

 //Function that calculates dt
double calculate_dt(CellTensor& U, const Grid& grid, const double& CFL) {

	Vector V;
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

void Calculate_Jacobians(const CellTensor& U, Vector& U_inlet, const BoundaryConditions& BoundaryTypes, const Grid& grid,
	iFaceTensor& i_Fluxes, jFaceTensor& j_Fluxes, iFaceTesseract& i_plus_Jacobians, iFaceTesseract& i_minus_Jacobians, jFaceTesseract& j_plus_Jacobians, jFaceTesseract& j_minus_Jacobians) {

	static Vector Ui, Uii, Up, Um, V, Vi, Vii, V1_Plus, V2_Plus, V1_Minus, V2_Minus, n, m;   
	double g = 5.72;
	double weight, dp, a, rho, u, v, p, nx, ny, uprime, pe, pp, h0, lp, lcp, ltp, lm, lcm, ltm, k; 
	pe = (gamma - 1);  
	Matrix X, Y, E;     

	// Calculate Jacobians and Explicit fluxes for i-faces on left boundary 
	for (int j = 0; j < Ny; ++j) { 
				
		Uii = U[0][j]; 
		tie(Ui, E) = Boundary2D(BoundaryTypes.left, Uii, U_inlet, grid.iNorms(0, j));    	 
		  
		Vi = constoPrim(Ui);   
		Vii = constoPrim(Uii);

		nx = grid.iNorms(0, j).x;   
		ny = grid.iNorms(0, j).y;   

		dp = fabs(Vii[3] - Vi[3]) / min(Vii[3], Vi[3]);  
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

		Up = weight * Ui + (1 - weight) * Uii;    
		Um = (1 - weight) * Ui + weight * Uii;  
		
		Vi = constoPrim(Up);
		Vii = constoPrim(Um);

		rho = Vi[0];  
		u = Vi[1];  
		v = Vi[2];
		p = Vi[3];		

		a = sqrt(gamma * p / rho);		 
		k = 1 / (a * a); 
		uprime = u * nx + v * ny;

		pp = 0.5 * (gamma - 1) * (u * u + v * v); 
		h0 = (Up[3] + p) / rho; 

		lp = 0.5 * (uprime + fabs(uprime));
		lcp = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) + 0.5 * (uprime - a + fabs(uprime - a)) - lp); 
		ltp = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) - 0.5 * (uprime - a + fabs(uprime - a))); 

		V1_Plus = { lcp * k,
				(u * lcp + a * nx * ltp) * k,
				(v * lcp + a * ny * ltp) * k,
				(h0 * lcp + a * uprime * ltp) * k };

		V2_Plus = { ltp / a, 
			u* ltp / a + nx * lcp,
			v* ltp / a + ny * lcp,
			h0* ltp / a + uprime * lcp };

		m = { pp, -u * pe, -v * pe, pe }; 
		n = { -uprime, nx, ny, 0 }; 		
		
		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(); 

		rho = Vii[0]; 
		u = Vii[1]; 
		v = Vii[2]; 
		p = Vii[3];  

		a = sqrt(gamma * p / rho); 
		k = 1 / (a * a); 
		uprime = u * nx + v * ny; 

		pp = 0.5 * (gamma - 1) * (u * u + v * v); 
		h0 = (Um[3] + p) / rho;  

		m = { pp, -u * pe, -v * pe, pe };
		n = { -uprime, nx, ny, 0 };

		lm = 0.5 * (uprime - fabs(uprime)); 
		lcm = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) + 0.5 * (uprime - a - fabs(uprime - a)) - lm); 
		ltm = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) - 0.5 * (uprime - a - fabs(uprime - a)));

		V1_Minus = { lcm * k,
					(u * lcm + a * nx * ltm) * k,
					(v * lcm + a * ny * ltm) * k,
					(h0 * lcm + a * uprime * ltm) * k };

		V2_Minus = { ltm / a, 
			u* ltm / a + nx * lcm,
			v* ltm / a + ny * lcm,
			h0* ltm / a + uprime * lcm }; 

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity();   
		i_plus_Jacobians[0][j] = X;  
		i_minus_Jacobians[0][j] = Y;  
		i_Fluxes[0][j] = X * Ui + Y * Uii;    
	}

	// Calculate Jacobians and Explicit fluxes for i-faces on right boundary
	for (int j = 0; j < Ny; ++j) {

		Ui = U[Nx - 1][j]; 
		tie(Uii, E) = Boundary2D(BoundaryTypes.right, Ui, U_inlet, grid.iNorms(Nx, j));     

		Vi = constoPrim(Ui);    
		Vii = constoPrim(Uii); 

		nx = grid.iNorms(Nx, j).x; 
		ny = grid.iNorms(Nx, j).y; 

		dp = fabs(Vii[3] - Vi[3]) / min(Vii[3], Vi[3]);
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

		Up = weight * Ui + (1 - weight) * Uii;
		Um = (1 - weight) * Ui + weight * Uii;

		Vi = constoPrim(Up);
		Vii = constoPrim(Um);

		rho = Vi[0];
		u = Vi[1];
		v = Vi[2];
		p = Vi[3];

		a = sqrt(gamma * p / rho);
		k = 1 / (a * a);
		uprime = u * nx + v * ny;

		pp = 0.5 * (gamma - 1) * (u * u + v * v);
		h0 = (Up[3] + p) / rho;

		lp = 0.5 * (uprime + fabs(uprime));
		lcp = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) + 0.5 * (uprime - a + fabs(uprime - a)) - lp);
		ltp = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) - 0.5 * (uprime - a + fabs(uprime - a)));

		V1_Plus = { lcp * k,
				(u * lcp + a * nx * ltp) * k,
				(v * lcp + a * ny * ltp) * k,
				(h0 * lcp + a * uprime * ltp) * k };

		V2_Plus = { ltp / a,
			u * ltp / a + nx * lcp,
			v * ltp / a + ny * lcp,
			h0 * ltp / a + uprime * lcp };

		m = { pp, -u * pe, -v * pe, pe };
		n = { -uprime, nx, ny, 0 };

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(); 

		rho = Vii[0];
		u = Vii[1];
		v = Vii[2];
		p = Vii[3];

		a = sqrt(gamma * p / rho);
		k = 1 / (a * a);
		uprime = u * nx + v * ny;

		pp = 0.5 * (gamma - 1) * (u * u + v * v);
		h0 = (Um[3] + p) / rho;

		m = { pp, -u * pe, -v * pe, pe };
		n = { -uprime, nx, ny, 0 };

		lm = 0.5 * (uprime - fabs(uprime));
		lcm = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) + 0.5 * (uprime - a - fabs(uprime - a)) - lm);
		ltm = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) - 0.5 * (uprime - a - fabs(uprime - a)));

		V1_Minus = { lcm * k,
				(u * lcm + a * nx * ltm) * k,
				(v * lcm + a * ny * ltm) * k,
				(h0 * lcm + a * uprime * ltm) * k };

		V2_Minus = { ltm / a,
			u * ltm / a + nx * lcm,
			v * ltm / a + ny * lcm,
			h0 * ltm / a + uprime * lcm };

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity();
		i_plus_Jacobians[Nx][j] = X; 
		i_minus_Jacobians[Nx][j] = Y; 
		i_Fluxes[Nx][j] = X * Ui + Y * Uii;   
	}

	// Calculate Jacobians and Explicit fluxes for j-faces on bottom boundary
	for (int i = 0; i < Nx; ++i) { 

		Uii = U[i][0];
		tie(Ui, E) = Boundary2D(BoundaryTypes.bottom, Uii, U_inlet, grid.jNorms(i, 0)); 

		Vi = constoPrim(Ui);  
		Vii = constoPrim(Uii);  

		nx = grid.jNorms(i, 0).x;  
		ny = grid.jNorms(i, 0).y;  

		dp = fabs(Vii[3] - Vi[3]) / min(Vii[3], Vi[3]);
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

		Up = weight * Ui + (1 - weight) * Uii;
		Um = (1 - weight) * Ui + weight * Uii;

		Vi = constoPrim(Up);
		Vii = constoPrim(Um);

		rho = Vi[0];
		u = Vi[1];
		v = Vi[2];
		p = Vi[3];

		a = sqrt(gamma * p / rho);
		k = 1 / (a * a);
		uprime = u * nx + v * ny;

		pp = 0.5 * (gamma - 1) * (u * u + v * v);
		h0 = (Up[3] + p) / rho;  

		lp = 0.5 * (uprime + fabs(uprime));
		lcp = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) + 0.5 * (uprime - a + fabs(uprime - a)) - lp);
		ltp = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) - 0.5 * (uprime - a + fabs(uprime - a)));

		V1_Plus = { lcp * k,
				(u * lcp + a * nx * ltp) * k,
				(v * lcp + a * ny * ltp) * k,
				(h0 * lcp + a * uprime * ltp) * k };

		V2_Plus = { ltp / a,
			u * ltp / a + nx * lcp,
			v * ltp / a + ny * lcp,
			h0 * ltp / a + uprime * lcp };

		m = { pp, -u * pe, -v * pe, pe };
		n = { -uprime, nx, ny, 0 };

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(); 

		rho = Vii[0];
		u = Vii[1];
		v = Vii[2];
		p = Vii[3];

		a = sqrt(gamma * p / rho);
		k = 1 / (a * a);
		uprime = u * nx + v * ny;

		pp = 0.5 * (gamma - 1) * (u * u + v * v);
		h0 = (Um[3] + p) / rho;

		m = { pp, -u * pe, -v * pe, pe };
		n = { -uprime, nx, ny, 0 };

		lm = 0.5 * (uprime - fabs(uprime));
		lcm = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) + 0.5 * (uprime - a - fabs(uprime - a)) - lm);
		ltm = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) - 0.5 * (uprime - a - fabs(uprime - a)));

		V1_Minus = { lcm * k,
				(u * lcm + a * nx * ltm) * k,
				(v * lcm + a * ny * ltm) * k,
				(h0 * lcm + a * uprime * ltm) * k };

		V2_Minus = { ltm / a,
			u * ltm / a + nx * lcm,
			v * ltm / a + ny * lcm,
			h0 * ltm / a + uprime * lcm };

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity();
		j_plus_Jacobians[i][0] = X;
		j_minus_Jacobians[i][0] = Y; 
  
		j_Fluxes[i][0] = X * Ui + Y * Uii;   

	}

	// Calculate Jacobians and Explicit fluxes for j-faces on top boundary
	for (int i = 0; i < Nx; ++i) {

		Ui = U[i][Ny - 1]; 
		tie(Uii, E) = Boundary2D(BoundaryTypes.top, Ui, U_inlet, grid.jNorms(i, Ny));  

		Vi = constoPrim(Ui); 
		Vii = constoPrim(Uii);  

		nx = grid.jNorms(i, Ny).x;   
		ny = grid.jNorms(i, Ny).y;   

		dp = fabs(Vii[3] - Vi[3]) / min(Vii[3], Vi[3]); 
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

		Up = weight * Ui + (1 - weight) * Uii; 
		Um = (1 - weight) * Ui + weight * Uii; 

		Vi = constoPrim(Up);
		Vii = constoPrim(Um);

		rho = Vi[0]; 
		u = Vi[1]; 
		v = Vi[2]; 
		p = Vi[3]; 

		a = sqrt(gamma * p / rho);
		k = 1 / (a * a);
		uprime = u * nx + v * ny;

		pp = 0.5 * (gamma - 1) * (u * u + v * v);
		h0 = (Up[3] + p) / rho;

		lp = 0.5 * (uprime + fabs(uprime));
		lcp = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) + 0.5 * (uprime - a + fabs(uprime - a)) - lp);
		ltp = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) - 0.5 * (uprime - a + fabs(uprime - a)));

		V1_Plus = { lcp * k,
				(u * lcp + a * nx * ltp) * k,
				(v * lcp + a * ny * ltp) * k,
				(h0 * lcp + a * uprime * ltp) * k };

		V2_Plus = { ltp / a,
			u * ltp / a + nx * lcp,
			v * ltp / a + ny * lcp,
			h0 * ltp / a + uprime * lcp };

		m = { pp, -u * pe, -v * pe, pe };
		n = { -uprime, nx, ny, 0 };

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(); 

		rho = Vii[0];
		u = Vii[1];
		v = Vii[2];
		p = Vii[3];

		a = sqrt(gamma * p / rho);
		k = 1 / (a * a);
		uprime = u * nx + v * ny;

		pp = 0.5 * (gamma - 1) * (u * u + v * v);
		h0 = (Um[3] + p) / rho;

		m = { pp, -u * pe, -v * pe, pe };
		n = { -uprime, nx, ny, 0 };

		lm = 0.5 * (uprime - fabs(uprime));
		lcm = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) + 0.5 * (uprime - a - fabs(uprime - a)) - lm);
		ltm = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) - 0.5 * (uprime - a - fabs(uprime - a)));

		V1_Minus = { lcm * k,
				(u * lcm + a * nx * ltm) * k,
				(v * lcm + a * ny * ltm) * k,
				(h0 * lcm + a * uprime * ltm) * k };

		V2_Minus = { ltm / a,
			u * ltm / a + nx * lcm,
			v * ltm / a + ny * lcm,
			h0 * ltm / a + uprime * lcm };

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(); 
		j_plus_Jacobians[i][Ny] = X; 
		j_minus_Jacobians[i][Ny] = Y; 
		j_Fluxes[i][Ny] = X * Ui + Y * Uii;   

	}

	// Inner i-faces  
	for (int i = 1; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {	

			Ui = U[i - 1][j];
			Uii = U[i][j]; 

			Vi = constoPrim(Ui);
			Vii = constoPrim(Uii);

			nx = grid.iNorms(i, j).x;  
			ny = grid.iNorms(i, j).y; 

			dp = fabs(Vii[3] - Vi[3]) / min(Vii[3], Vi[3]); 
			weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

			Up = weight * Ui + (1 - weight) * Uii; 
			Um = (1 - weight) * Ui + weight * Uii; 

			Vi = constoPrim(Up);
			Vii = constoPrim(Um);

			rho = Vi[0]; 
			u = Vi[1]; 
			v = Vi[2]; 
			p = Vi[3]; 

			a = sqrt(gamma * p / rho);
			k = 1 / (a * a);  
			uprime = u * nx + v * ny;

			pp = 0.5 * (gamma - 1) * (u * u + v * v);
			h0 = (Up[3] + p) / rho;

			lp = 0.5 * (uprime + fabs(uprime));
			lcp = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) + 0.5 * (uprime - a + fabs(uprime - a)) - lp);
			ltp = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) - 0.5 * (uprime - a + fabs(uprime - a)));

			V1_Plus = { lcp * k, 
				(u * lcp + a * nx * ltp) * k, 
				(v * lcp + a * ny * ltp)* k, 
				(h0 * lcp + a * uprime * ltp)* k }; 

			V2_Plus = { ltp / a,
				u * ltp / a + nx * lcp,
				v * ltp / a + ny * lcp,
				h0 * ltp / a + uprime * lcp };

			m = { pp, -u * pe, -v * pe, pe };
			n = { -uprime, nx, ny, 0 };

			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity();
			 
			rho = Vii[0];
			u = Vii[1];
			v = Vii[2];
			p = Vii[3];

			a = sqrt(gamma * p / rho);
			k = 1 / (a * a);
			uprime = u * nx + v * ny;

			pp = 0.5 * (gamma - 1) * (u * u + v * v);
			h0 = (Um[3] + p) / rho;

			m = { pp, -u * pe, -v * pe, pe };
			n = { -uprime, nx, ny, 0 };

			lm = 0.5 * (uprime - fabs(uprime));
			lcm = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) + 0.5 * (uprime - a - fabs(uprime - a)) - lm);
			ltm = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) - 0.5 * (uprime - a - fabs(uprime - a)));

			V1_Minus = { lcm * k, 
				(u * lcm + a * nx * ltm) * k,
				(v * lcm + a * ny * ltm) * k,
				(h0 * lcm + a * uprime * ltm) * k };

			V2_Minus = { ltm / a,
				u * ltm / a + nx * lcm,
				v * ltm / a + ny * lcm,
				h0 * ltm / a + uprime * lcm };
			
			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(); 
			 
			i_plus_Jacobians[i][j] = X;
			i_minus_Jacobians[i][j] = Y;
			i_Fluxes[i][j] = X * Ui + Y * Uii;   
		}
	}

	// Inner j-faces 
	for (int i = 0; i < Nx; ++i) {
		for (int j = 1; j < Ny; ++j) {
			
			Ui = U[i][j - 1]; 
			Uii = U[i][j];  
			Vi = constoPrim(Ui); 
			Vii = constoPrim(Uii);  

			nx = grid.jNorms(i, j).x;  
			ny = grid.jNorms(i, j).y;  

			dp = fabs(Vii[3] - Vi[3]) / min(Vii[3], Vi[3]); 
			weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

			Up = weight * Ui + (1 - weight) * Uii; 
			Um = (1 - weight) * Ui + weight * Uii; 

			Vi = constoPrim(Up);  
			Vii = constoPrim(Um);  

			rho = Vi[0]; 
			u = Vi[1]; 
			v = Vi[2]; 
			p = Vi[3]; 

			a = sqrt(gamma * p / rho);
			k = 1 / (a * a); 
			uprime = u * nx + v * ny;

			pp = 0.5 * (gamma - 1) * (u * u + v * v);
			h0 = (Up[3] + p) / rho;

			lp = 0.5 * (uprime + fabs(uprime));
			lcp = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) + 0.5 * (uprime - a + fabs(uprime - a)) - lp);
			ltp = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) - 0.5 * (uprime - a + fabs(uprime - a)));

			V1_Plus = { lcp * k,
				(u * lcp + a * nx * ltp)* k,
				(v * lcp + a * ny * ltp)* k,
				(h0 * lcp + a * uprime * ltp)* k };

			V2_Plus = { ltp / a,
				u * ltp / a + nx * lcp,
				v * ltp / a + ny * lcp,
				h0 * ltp / a + uprime * lcp };

			m = { pp, -u * pe, -v * pe, pe };
			n = { -uprime, nx, ny, 0 };

			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity();   

			rho = Vii[0];
			u = Vii[1];
			v = Vii[2];
			p = Vii[3];

			a = sqrt(gamma * p / rho);
			k = 1 / (a * a);
			uprime = u * nx + v * ny;

			pp = 0.5 * (gamma - 1) * (u * u + v * v);
			h0 = (Um[3] + p) / rho;

			m = { pp, -u * pe, -v * pe, pe };
			n = { -uprime, nx, ny, 0 };

			lm = 0.5 * (uprime - fabs(uprime));
			lcm = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) + 0.5 * (uprime - a - fabs(uprime - a)) - lm);
			ltm = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) - 0.5 * (uprime - a - fabs(uprime - a)));

			V1_Minus = { lcm * k, 
				(u * lcm + a * nx * ltm)* k, 
				(v * lcm + a * ny * ltm)* k, 
				(h0 * lcm + a * uprime * ltm)* k };  

			V2_Minus = { ltm / a,
				u * ltm / a + nx * lcm,
				v * ltm / a + ny * lcm,
				h0 * ltm / a + uprime * lcm };


			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(); 
			j_plus_Jacobians[i][j] = X; 
			j_minus_Jacobians[i][j] = Y;
			j_Fluxes[i][j] = X * Ui + Y * Uii;   
		}
	}

} 

void solveLeftLine(CellTensor& U, CellTensor& dU_new, Vector U_inlet, CellTensor& dU_old, const Grid& grid,
	BoundaryConditions BoundaryType, const int i, const double dt,
	iFaceTensor& i_Fluxes, jFaceTensor& j_Fluxes, iFaceTesseract& i_plus_Jacobians, iFaceTesseract& i_minus_Jacobians, jFaceTesseract& j_plus_Jacobians, jFaceTesseract& j_minus_Jacobians) {

	static Matrix A; // Relates to U(j+1) 
	static Matrix B; // Relates to U(j) 
	static Matrix C; // Relates to U(j-1) 
	static Vector F; // Right hand side  

	// Intermediate matrices
	static Matrix alpha;
	static array<array<double, 4>, Ny> v;
	static array<array<array<double,4>, 4>, Ny> g; 

	// Grab boundary values
	Duo bottomBC = Boundary2D(BoundaryType.bottom, U[i][0], U_inlet, grid.jNorms(i, 0));
	Duo topBC = Boundary2D(BoundaryType.top, U[i][Ny - 1], U_inlet, grid.jNorms(i, Ny));

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity()
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

		A = grid.Volume(i, j) / dt * identity()
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

	A = grid.Volume(i, 0) / dt * identity()
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

void solveMiddleLine(CellTensor& U, CellTensor& dU_new, Vector U_inlet, CellTensor& dU_old, const Grid& grid,   
	BoundaryConditions BoundaryType, const int i, const double dt,
	iFaceTensor& i_Fluxes, jFaceTensor& j_Fluxes, iFaceTesseract& i_plus_Jacobians, iFaceTesseract& i_minus_Jacobians, jFaceTesseract& j_plus_Jacobians, jFaceTesseract& j_minus_Jacobians) {

	//auto start = TIME;

	static Matrix A; // Relates to U(j+1) 
	static Matrix B; // Relates to U(j) 
	static Matrix C; // Relates to U(j-1) 
	static Vector F; // Right hand side  

	// Intermediate matrices
	static Matrix alpha; 
	static array<array<double, 4>, Ny> v;
	static array<array<array<double, 4>, 4>, Ny> g; 

	// Grab boundary values
	Duo bottomBC = Boundary2D(BoundaryType.bottom, U[i][0], U_inlet, grid.jNorms(i, 0)); 
	Duo topBC = Boundary2D(BoundaryType.top, U[i][Ny - 1], U_inlet, grid.jNorms(i, Ny));   

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity() 
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

		A = grid.Volume(i, j) / dt * identity() 
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

	A = grid.Volume(i, 0) / dt * identity()
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

	//auto end = TIME;
	//DURATION duration = end - start; 
	//cout << duration.count() << endl; 
}

void solveRightLine(CellTensor& U, CellTensor& dU_new, Vector U_inlet, CellTensor& dU_old, const Grid& grid,
	BoundaryConditions BoundaryType, const int i, const double dt,
	iFaceTensor& i_Fluxes, jFaceTensor& j_Fluxes, iFaceTesseract& i_plus_Jacobians, iFaceTesseract& i_minus_Jacobians, jFaceTesseract& j_plus_Jacobians, jFaceTesseract& j_minus_Jacobians) {

	static Matrix A; // Relates to U(j+1) 
	static Matrix B; // Relates to U(j) 
	static Matrix C; // Relates to U(j-1) 
	static Vector F; // Right hand side  

	// Intermediate matrices
	static Matrix alpha;
	static array<array<double, 4>, Ny> v; 
	static array<array<array<double, 4>, 4>, Ny> g; 

	// Grab boundary values
	Duo bottomBC = Boundary2D(BoundaryType.bottom, U[i][0], U_inlet, grid.jNorms(i, 0));
	Duo topBC = Boundary2D(BoundaryType.top, U[i][Ny - 1], U_inlet, grid.jNorms(i, Ny));

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity()
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

		A = grid.Volume(i, j) / dt * identity()
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

	A = grid.Volume(i, 0) / dt * identity()
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

double calculateInnerResidual(CellTensor& U, CellTensor& dU_new, CellTensor& dU_old, const Grid& grid, const double dt,
	iFaceTensor& i_Fluxes, jFaceTensor& j_Fluxes, iFaceTesseract& i_plus_Jacobians, iFaceTesseract& i_minus_Jacobians, jFaceTesseract& j_plus_Jacobians, jFaceTesseract& j_minus_Jacobians) {
	//
	double inner_residual = 0.0;
	double res = 0.0; 

	static Vector A; // Relates to U(j+1) 
	static Vector B; // Relates to U(j) 
	static Vector C; // Relates to U(j-1) 
	static double F; // Right hand side  

	Vector one = { 1, 0, 0, 0 };

	// Middle Boundary
	for (int i = 1; i < Nx - 1; ++i) {
		for (int j = 1; j < Ny - 1; ++j) {

			B = j_minus_Jacobians[i][j + 1][0] * (grid.jArea(i, j + 1)); 

			A = grid.Volume(i, j) / dt * one 
				- i_minus_Jacobians[i][j][0] * grid.iArea(i, j) 
				+ i_plus_Jacobians[i + 1][j][0] * grid.iArea(i + 1, j)
				- j_minus_Jacobians[i][j][0] * grid.jArea(i, j) 
				+ j_plus_Jacobians[i][j + 1][0] * grid.jArea(i, j + 1);

			C = j_plus_Jacobians[i][j][0] * (-grid.jArea(i, j)); 

			F = i_Fluxes[i][j][0] * (grid.iArea(i, j)) 
				+ i_Fluxes[i + 1][j][0] * (-grid.iArea(i + 1, j)) 
				+ j_Fluxes[i][j][0] * (grid.jArea(i, j)) 
				+ j_Fluxes[i][j + 1][0] * (-grid.jArea(i, j + 1)) 
				+ i_plus_Jacobians[i][j][0] * (grid.iArea(i, j)) * dU_old[i - 1][j]
				+ i_minus_Jacobians[i + 1][j][0] * (-grid.iArea(i + 1, j)) * dU_old[i + 1][j];   

			res = B * dU_new[i][j + 1] + A * dU_new[i][j] + C * dU_new[i][j - 1] - F;
			inner_residual += res * res;  

		}
	}
	return sqrt(inner_residual); 
}

double calculateResidual(const CellTensor& U, const Grid& grid, const iFaceTensor& i_Fluxes, const jFaceTensor& j_Fluxes) {

	double res = 0.0;
	double intres;

	for (int i = 1; i < Nx - 1; ++i) {
		for (int j = 1; j < Ny - 1; ++j) {
			intres = (i_Fluxes[i][j][0] * (-grid.iArea(i, j))			// Left Side		   
				+ i_Fluxes[i + 1][j][0] * (grid.iArea(i + 1, j))		// Right Side 
				+ j_Fluxes[i][j][0] * (-grid.jArea(i, j))			// Bottom Side  
				+ j_Fluxes[i][j + 1][0] * (grid.jArea(i, j + 1)))/ grid.Volume(i, j); 

			res = res + intres * intres;
		
		}
	}

	return sqrt(res); 
}