// 2D FVS Library.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "framework.h"
#include "2DFVSLibrary.h"


Solver::Solver(const inlet_conditions& INLET, Grid& grid, BoundaryConditions BoundaryType, const double CFL, const double Tw, const int& progress_update) : Tw(Tw), INLET(INLET), V_inlet(V_inlet), U_inlet(U_inlet), grid(grid), BoundaryType(BoundaryType), U(U), dU_new(dU_new), dU_old(dU_old), 
i_Fluxes(i_Fluxes), j_Fluxes(j_Fluxes), i_plus_inviscid_Jacobians(i_plus_inviscid_Jacobians), j_plus_inviscid_Jacobians(j_plus_inviscid_Jacobians), i_minus_inviscid_Jacobians(i_minus_inviscid_Jacobians), 
j_minus_inviscid_Jacobians(j_minus_inviscid_Jacobians), i_viscous_Jacobians(i_viscous_Jacobians), j_viscous_Jacobians(j_viscous_Jacobians),  gridtype(gridtype), dt(dt), CFL(CFL), inner_residual(inner_residual), outer_residual(outer_residual),
progress_update(progress_update){

	outer_residual = 1.0;
	inner_residual = 1.0; 

	V_inlet = { INLET.rho, INLET.u, INLET.v, INLET.p };

	gridtype;   
	if (dynamic_cast<RampGrid*>(&grid)) gridtype = "Ramp";  
	else if (dynamic_cast<CylinderGrid*>(&grid)) gridtype = "Cylinder";  
	else if (dynamic_cast<FlatPlateGrid*>(&grid)) gridtype = "Flat Plate";  
	else gridtype = "Unknown";  


	U_inlet = primtoCons(V_inlet); 

	for (int i = 0; i < Nx; ++i) { 
		for (int j = 0; j < Ny; ++j) { 

			U[i][j] = U_inlet; 

			for (int k = 0; k < 4; ++k) { 

				dU_new[i][j][k] = 0.0; 
				dU_old[i][j][k] = 0.0; 
			} 
		}
	}
};   

Vector Solver::primtoCons(const Vector& V) { 

	Vector U = zerosV(); 
	U[0] = V[0];
	U[1] = V[0] * V[1];
	U[2] = V[0] * V[2];
	U[3] = V[3] / (gamma - 1) + 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2]);

	return U;
}

Vector Solver::constoPrim(const Vector& U) { 

	static Vector V = zerosV(); 
	V[0] = U[0];
	V[1] = U[1] / U[0];
	V[2] = U[2] / U[0];
	V[3] = (U[3] - 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2])) * (gamma - 1);

	return V;
}

Vector Solver::constoViscPrim(const Vector& U) { 
	Vector result;
	result[0] = U[0];
	result[1] = U[1] / U[0];
	result[2] = U[2] / U[0];
	result[3] = computeTemperature(U); 
	return result;  
}

Matrix Solver::inviscid_boundary_2D_E(BoundaryCondition type, const Vector& U, const Point& normals) {     

	Matrix E = zerosM();  	

	switch (type) {

	case BoundaryCondition::Inlet: 
		return E;  

	case BoundaryCondition::Outlet: 
		
		return onesM();  

	case BoundaryCondition::IsothermalWall:  
	 
		E = { {{1, 0, 0, 0}, 
			  {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
			  { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
			  { 0, 0, 0, 1 }}
		};

		return E;

	case BoundaryCondition::AdiabaticWall:

		E = { {{1, 0, 0, 0}, 
			  {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
			  { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
			  { 0, 0, 0, 1 }}
		};

		return E;

	case BoundaryCondition::Symmetry: 		

		E = { {{1, 0, 0, 0},
			  {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
			  { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
			  { 0, 0, 0, 1 }}
		};

		return E; 
	default:
		throw invalid_argument("Unknown boundary condition type.");
	}

}

Vector Solver::inviscid_boundary_2D_U(BoundaryCondition type, const Vector& U, const Point& normals) {

	Vector ghost;  

	double u = 0.0, v = 0.0;

	switch (type) {

	case BoundaryCondition::Inlet:
		return U_inlet; 

	case BoundaryCondition::Outlet:
		
		return U;

	case BoundaryCondition::IsothermalWall: 

		u = U[1] / U[0];
		v = U[2] / U[0];

		ghost[0] = U[0];
		ghost[1] = U[0] * (u - 2 * (u * normals.x + v * normals.y) * normals.x);
		ghost[2] = U[0] * (v - 2 * (u * normals.x + v * normals.y) * normals.y);
		ghost[3] = U[3];

		return ghost;

	case BoundaryCondition::AdiabaticWall:

		u = U[1] / U[0];
		v = U[2] / U[0];

		ghost[0] = U[0];
		ghost[1] = U[0] * (u - 2 * (u * normals.x + v * normals.y) * normals.x);
		ghost[2] = U[0] * (v - 2 * (u * normals.x + v * normals.y) * normals.y);
		ghost[3] = U[3];

		return ghost;

	case BoundaryCondition::Symmetry:

		u = U[1] / U[0]; 
		v = U[2] / U[0]; 

		ghost[0] = U[0]; 
		ghost[1] = U[0] * (u - 2 * (u * normals.x + v * normals.y) * normals.x);
		ghost[2] = U[0] * (v - 2 * (u * normals.x + v * normals.y) * normals.y);
		ghost[3] = U[3];

		return ghost;

	default:
		throw invalid_argument("Unknown boundary condition type.");
	}

}

Matrix Solver::viscous_boundary_2D_E(BoundaryCondition type, const Vector& U, const Point& normals) { 

	Matrix E = zerosM();
	double Ti; 

	switch (type) {

	case BoundaryCondition::Inlet:
		return E; 

	case BoundaryCondition::Outlet:

		return onesM();

	case BoundaryCondition::IsothermalWall:

		Ti = computeTemperature(U); 

		E = { {{2 * Ti/Tw - 1 , 0, 0, 2 * U[0]/Tw}, 
			  { 0, -1, 0, 0 },
			  { 0, 0, -1, 0},
			  { 0, 0, 0, -1 }}
		};

		return E;

	case BoundaryCondition::AdiabaticWall:

		Ti = computeTemperature(U);

		E = { {{2 * Ti / Tw - 1 , 0, 0, 2 * U[0] / Tw}, 
			  { 0, -1, 0, 0 },
			  { 0, 0, -1, 0},
			  { 0, 0, 0, 1 }}
		};

		return E; 

	case BoundaryCondition::Symmetry:

		E = { {{1, 0, 0, 0},
			  {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
			  { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
			  { 0, 0, 0, 1 }}
		};

		return E;
	default:
		throw invalid_argument("Unknown boundary condition type.");
	}

}

Vector Solver::viscous_boundary_2D_U(BoundaryCondition type, const Vector& U, const Point& normals) {

	Vector ghost;
	double T_in; 
	double u = 0.0, v = 0.0; 

	switch (type) {

	case BoundaryCondition::Inlet:
		return U_inlet;

	case BoundaryCondition::Outlet:

		return U;

	case BoundaryCondition::IsothermalWall:

		T_in = computeTemperature(U);
		
		ghost[0] = 2 * (U[0] * T_in / Tw) - U[0]; 
		ghost[1] = -U[1];  
		ghost[2] = -U[2]; 
		ghost[3] = U[3];

		return ghost;

	case BoundaryCondition::AdiabaticWall:

		T_in = Tw;   

		ghost[0] = 2 * (U[0] * T_in / Tw) - U[0]; 
		ghost[1] = -U[1];  
		ghost[2] = -U[2];  
		ghost[3] = U[3]; 

		return ghost;

	case BoundaryCondition::Symmetry:

		u = U[1] / U[0]; 
		v = U[2] / U[0]; 

		ghost[0] = U[0];
		ghost[1] = U[0] * (u - 2 * (u * normals.x + v * normals.y) * normals.x);
		ghost[2] = U[0] * (v - 2 * (u * normals.x + v * normals.y) * normals.y);
		ghost[3] = U[3];

		return ghost;

	default:
		throw invalid_argument("Unknown boundary condition type.");
	}

}

void Solver::compute_dt() {  

	Vector V;
	double dx, dy, c, dt_old; 
	dt = 1; 

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
}

void Solver::solve_inviscid () { 

	cout << "Running Inviscid DPLR for " << Nx << " by " << Ny << " " << gridtype << " with a CFL of " << fixed << setprecision(2) << CFL << "..." << "\n\n";
	string filename = "Inviscid " + to_string(Nx) + "x" + to_string(Ny) + "_" + gridtype + "_Solution.csv"; 
	
	auto start = TIME; 
	int counter = 0;

	while (outer_residual >= 1e-6) { 

		compute_dt();

		solve_inviscid_timestep();

		compute_outer_residual();

		if (counter == 0) outer_residual = 1.0;
		counter++;

		if (counter % progress_update == 0) { 
			auto end = TIME;
			DURATION duration = end - start;
			cout << "Iteration: " << counter << "\t Residual: " << fixed << scientific << setprecision(3) << outer_residual 
				 << fixed << setprecision(2) << "\tElapsed time: " << duration.count() << " seconds" << endl;
		}

		if (counter % 1000 == 0) { 
			write_2d_csv(filename); 
		}


	}

	write_2d_csv(filename);  

	auto end = TIME;
	DURATION duration = end - start; 
	cout << "Program complete in " << duration.count() << endl; 


}  

void Solver::solve_viscous() {

	cout << "Running Viscous DPLR for " << Nx << " by " << Ny << " " << gridtype << " with a CFL of " << fixed << setprecision(2) << CFL << "..." << "\n\n"; 
	string filename = "Viscous " + to_string(Nx) + "x" + to_string(Ny) + "_" + gridtype + "_Solution.csv";

	auto start = TIME;
	int counter = 0;

	while (outer_residual >= 1e-6) {

		compute_dt();

		solve_viscous_timestep();  

		compute_outer_residual();

		if (counter == 0) outer_residual = 1.0;
		counter++;

		if (counter % progress_update == 0) { 
			auto end = TIME;
			DURATION duration = end - start;
			cout << "Iteration: " << counter << "\t Residual: " << fixed << scientific << setprecision(3) << outer_residual
				<< fixed << setprecision(2) << "\tElapsed time: " << duration.count() << " seconds" << endl;
		}

		if (counter % 1000 == 0) {
			write_2d_csv(filename);
		}

	}

	write_2d_csv(filename); 
	write_1d_csv("1D VISCOUS.csv");     

	auto end = TIME;
	DURATION duration = end - start;
	cout << "Program complete in " << duration.count() << endl;



}

void Solver::solve_viscous_timestep() {

	inner_residual = 1.0;

	while (inner_residual >= 1e-8) {

		compute_viscous_jacobians();  

		solve_left_line_viscous();

		for (int i = 1; i < Nx - 1; ++i) {
			solve_middle_line_viscous(i); 
		}

		solve_right_line_viscous(); 

		compute_inner_residual();
		swap(dU_old, dU_new);
	}

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < 4; ++k) {
				U[i][j][k] += dU_old[i][j][k];
			}
		}
	}

}

void Solver::solve_inviscid_timestep() {

	inner_residual = 1.0; 

	while (inner_residual >= 1e-8) {	

		compute_inviscid_jacobians(); 

		solve_left_line_inviscid();

		for (int i = 1; i < Nx - 1; ++i) {
			solve_middle_line_inviscid(i);
		}

		solve_right_line_inviscid();
			
		compute_inner_residual(); 
		swap(dU_old, dU_new); 
	}

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < 4; ++k) {
				U[i][j][k] += dU_old[i][j][k];    
			}
		}
	}

}

void Solver::compute_inviscid_jacobians() {

	static Vector Ui, Uii, Up, Um, V1_Plus, V2_Plus, V1_Minus, V2_Minus, n, m;   
	Inviscid_State Si, Sii; 
	double g = 5.72;
	double weight, dp, pi, pii, pe, nx, ny, l, lc, lt;
	pe = (gamma - 1);  
	Matrix X, Y;     

	// Calculate Jacobians and Explicit fluxes for i-faces on left boundary 
	for (int j = 0; j < Ny; ++j) { 
				
		Uii = U[0][j]; 
		Ui = inviscid_boundary_2D_U(BoundaryType.left, Uii, grid.iNorms(0, j));    	 

		nx = grid.iNorms(0, j).x;   
		ny = grid.iNorms(0, j).y;   

		pi = computePressure(Ui);
		pii = computePressure(Uii);
		dp = fabs(pii - pi) / min(pii, pi);
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

		Up = weight * Ui + (1 - weight) * Uii;
		Um = (1 - weight) * Ui + weight * Uii;
		
		Si = compute_inviscid_state(Up, nx, ny);
		Sii = compute_inviscid_state(Um, nx, ny);

		l = 0.5 * (Si.uprime + fabs(Si.uprime));
		lc = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - l);
		lt = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));

		V1_Plus = { lc * Si.k,
			(Si.u * lc + Si.a * nx * lt) * Si.k,
			(Si.v * lc + Si.a * ny * lt) * Si.k,
			(Si.h0 * lc + Si.a * Si.uprime * lt) * Si.k };

		V2_Plus = { lt / Si.a,
			Si.u * lt / Si.a + nx * lc,
			Si.v * lt / Si.a + ny * lc,
			Si.h0 * lt / Si.a + Si.uprime * lc };

		m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
		n = { -Si.uprime, nx, ny, 0 };

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity();

		m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
		n = { -Sii.uprime, nx, ny, 0 };

		l = 0.5 * (Sii.uprime - fabs(Sii.uprime));
		lc = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - l);
		lt = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

		V1_Minus = { lc * Sii.k,
				(Sii.u * lc + Sii.a * nx * lt) * Sii.k,
				(Sii.v * lc + Sii.a * ny * lt) * Sii.k,
				(Sii.h0 * lc + Sii.a * Sii.uprime * lt) * Sii.k };

		V2_Minus = { lt / Sii.a,
			Sii.u * lt / Sii.a + nx * lc,
			Sii.v * lt / Sii.a + ny * lc,
			Sii.h0 * lt / Sii.a + Sii.uprime * lc };


		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity();   
		i_plus_inviscid_Jacobians[0][j] = X;  
		i_minus_inviscid_Jacobians[0][j] = Y;  
		i_Fluxes[0][j] = X * Ui + Y * Uii;    
	}




	// Calculate Jacobians and Explicit fluxes for i-faces on right boundary
	for (int j = 0; j < Ny; ++j) {

		Ui = U[Nx - 1][j]; 
		Uii = inviscid_boundary_2D_U(BoundaryType.right, Ui, grid.iNorms(Nx, j));     

		nx = grid.iNorms(Nx, j).x; 
		ny = grid.iNorms(Nx, j).y; 

		pi = computePressure(Ui);
		pii = computePressure(Uii);
		dp = fabs(pii - pi) / min(pii, pi);
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

		Up = weight * Ui + (1 - weight) * Uii;
		Um = (1 - weight) * Ui + weight * Uii;

		Si = compute_inviscid_state(Up, nx, ny);
		Sii = compute_inviscid_state(Um, nx, ny);

		l = 0.5 * (Si.uprime + fabs(Si.uprime));
		lc = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - l);
		lt = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));

		V1_Plus = { lc * Si.k,
			(Si.u * lc + Si.a * nx * lt) * Si.k,
			(Si.v * lc + Si.a * ny * lt) * Si.k,
			(Si.h0 * lc + Si.a * Si.uprime * lt) * Si.k };

		V2_Plus = { lt / Si.a,
			Si.u * lt / Si.a + nx * lc,
			Si.v * lt / Si.a + ny * lc,
			Si.h0 * lt / Si.a + Si.uprime * lc };

		m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
		n = { -Si.uprime, nx, ny, 0 };

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity();

		m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
		n = { -Sii.uprime, nx, ny, 0 };

		l = 0.5 * (Sii.uprime - fabs(Sii.uprime));
		lc = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - l);
		lt = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

		V1_Minus = { lc * Sii.k,
				(Sii.u * lc + Sii.a * nx * lt) * Sii.k,
				(Sii.v * lc + Sii.a * ny * lt) * Sii.k,
				(Sii.h0 * lc + Sii.a * Sii.uprime * lt) * Sii.k };

		V2_Minus = { lt / Sii.a,
			Sii.u * lt / Sii.a + nx * lc,
			Sii.v * lt / Sii.a + ny * lc,
			Sii.h0 * lt / Sii.a + Sii.uprime * lc };

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity();
		

		i_plus_inviscid_Jacobians[Nx][j] = X; 
		i_minus_inviscid_Jacobians[Nx][j] = Y; 
		i_Fluxes[Nx][j] = X * Ui + Y * Uii;   
	}

	// Calculate Jacobians and Explicit fluxes for j-faces on bottom boundary
	for (int i = 0; i < Nx; ++i) { 

		Uii = U[i][0];
		Ui = inviscid_boundary_2D_U(BoundaryType.bottom, Uii, grid.jNorms(i, 0)); 

		nx = grid.jNorms(i, 0).x;  
		ny = grid.jNorms(i, 0).y;  

		pi = computePressure(Ui); 
		pii = computePressure(Uii); 
		dp = fabs(pii - pi) / min(pii, pi); 
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1)); 


		Up = weight * Ui + (1 - weight) * Uii; 
		Um = (1 - weight) * Ui + weight * Uii; 

		Si = compute_inviscid_state(Up, nx, ny);
		Sii = compute_inviscid_state(Um, nx, ny);

		l = 0.5 * (Si.uprime + fabs(Si.uprime));
		lc = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - l);
		lt = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));

		V1_Plus = { lc * Si.k,
				(Si.u * lc + Si.a * nx * lt) * Si.k,
				(Si.v * lc + Si.a * ny * lt) * Si.k,
				(Si.h0 * lc + Si.a * Si.uprime * lt) * Si.k };

		V2_Plus = { lt / Si.a,
			Si.u * lt / Si.a + nx * lc,
			Si.v * lt / Si.a + ny * lc,
			Si.h0 * lt / Si.a + Si.uprime * lc };

		m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
		n = { -Si.uprime, nx, ny, 0 };

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity();

		m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
		n = { -Sii.uprime, nx, ny, 0 };

		l = 0.5 * (Sii.uprime - fabs(Sii.uprime));
		lc = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - l);
		lt = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

		V1_Minus = { lc * Sii.k,
				(Sii.u * lc + Sii.a * nx * lt) * Sii.k,
				(Sii.v * lc + Sii.a * ny * lt) * Sii.k,
				(Sii.h0 * lc + Sii.a * Sii.uprime * lt) * Sii.k };

		V2_Minus = { lt / Sii.a,
			Sii.u * lt / Sii.a + nx * lc,
			Sii.v * lt / Sii.a + ny * lc,
			Sii.h0 * lt / Sii.a + Sii.uprime * lc };

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity();

		j_plus_inviscid_Jacobians[i][0] = X;
		j_minus_inviscid_Jacobians[i][0] = Y;   
		j_Fluxes[i][0] = X * Ui + Y * Uii;   

	}

	// Calculate Jacobians and Explicit fluxes for j-faces on top boundary
	for (int i = 0; i < Nx; ++i) {

		Ui = U[i][Ny - 1]; 
		Uii = inviscid_boundary_2D_U(BoundaryType.top, Ui, grid.jNorms(i, Ny));  

		nx = grid.jNorms(i, Ny).x;   
		ny = grid.jNorms(i, Ny).y;   

		pi = computePressure(Ui); 
		pii = computePressure(Uii); 
		dp = fabs(pii - pi) / min(pii, pi); 
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

		Up = weight * Ui + (1 - weight) * Uii;
		Um = (1 - weight) * Ui + weight * Uii;

		Si = compute_inviscid_state(Up, nx, ny);
		Sii = compute_inviscid_state(Um, nx, ny);

		l = 0.5 * (Si.uprime + fabs(Si.uprime));
		lc = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - l);
		lt = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));

		V1_Plus = { lc * Si.k,
				(Si.u * lc + Si.a * nx * lt) * Si.k,
				(Si.v * lc + Si.a * ny * lt) * Si.k,
				(Si.h0 * lc + Si.a * Si.uprime * lt) * Si.k };

		V2_Plus = { lt / Si.a,
			Si.u * lt / Si.a + nx * lc,
			Si.v * lt / Si.a + ny * lc,
			Si.h0 * lt / Si.a + Si.uprime * lc };

		m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
		n = { -Si.uprime, nx, ny, 0 };

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity();

		m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
		n = { -Sii.uprime, nx, ny, 0 };

		l = 0.5 * (Sii.uprime - fabs(Sii.uprime));
		lc = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - l);
		lt = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

		V1_Minus = { lc * Sii.k,
				(Sii.u * lc + Sii.a * nx * lt) * Sii.k,
				(Sii.v * lc + Sii.a * ny * lt) * Sii.k,
				(Sii.h0 * lc + Sii.a * Sii.uprime * lt) * Sii.k };

		V2_Minus = { lt / Sii.a,
			Sii.u * lt / Sii.a + nx * lc,
			Sii.v * lt / Sii.a + ny * lc,
			Sii.h0 * lt / Sii.a + Sii.uprime * lc };

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity(); 

		j_plus_inviscid_Jacobians[i][Ny] = X; 
		j_minus_inviscid_Jacobians[i][Ny] = Y; 
		j_Fluxes[i][Ny] = X * Ui + Y * Uii;   

	}

	// Inner i-faces  
	for (int i = 1; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {	

			nx = grid.iNorms(i, j).x;
			ny = grid.iNorms(i, j).y;

			Ui = U[i - 1][j];
			Uii = U[i][j]; 

			pi = computePressure(Ui);
			pii = computePressure(Uii); 
			dp = fabs(pii - pi) / min(pii, pi);
			weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
							

			Up = weight * Ui + (1 - weight) * Uii; 
			Um = (1 - weight) * Ui + weight * Uii; 

			Si = compute_inviscid_state(Up, nx, ny); 
			Sii = compute_inviscid_state(Um, nx, ny);

			l = 0.5 * (Si.uprime + fabs(Si.uprime));
			lc = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - l);
			lt = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));

			V1_Plus = { lc * Si.k,
				(Si.u * lc + Si.a * nx * lt) * Si.k,
				(Si.v * lc + Si.a * ny * lt) * Si.k,
				(Si.h0 * lc + Si.a * Si.uprime * lt) * Si.k };

			V2_Plus = { lt / Si.a,
				Si.u * lt / Si.a + nx * lc,
				Si.v * lt / Si.a + ny * lc,
				Si.h0 * lt / Si.a + Si.uprime * lc };

			m = { Si.pp, -Si.u * pe, -Si.v * pe, pe }; 
			n = { -Si.uprime, nx, ny, 0 }; 

			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity();
			

			m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe }; 
			n = { -Sii.uprime, nx, ny, 0 }; 

			l = 0.5 * (Sii.uprime - fabs(Sii.uprime));
			lc = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - l);
			lt = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

			V1_Minus = { lc * Sii.k,
				(Sii.u * lc + Sii.a * nx * lt) * Sii.k,
				(Sii.v * lc + Sii.a * ny * lt) * Sii.k,
				(Sii.h0 * lc + Sii.a * Sii.uprime * lt) * Sii.k };

			V2_Minus = { lt / Sii.a,
				Sii.u * lt / Sii.a + nx * lc,
				Sii.v * lt / Sii.a + ny * lc,
				Sii.h0 * lt / Sii.a + Sii.uprime * lc };
			
			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity(); 
			 
			i_plus_inviscid_Jacobians[i][j] = X;
			i_minus_inviscid_Jacobians[i][j] = Y;
			i_Fluxes[i][j] = X * Ui + Y * Uii;   
		}
	}




	// Inner j-faces 
	for (int i = 0; i < Nx; ++i) {
		for (int j = 1; j < Ny; ++j) {
			
			nx = grid.jNorms(i, j).x;
			ny = grid.jNorms(i, j).y;

			Ui = U[i][j - 1]; 
			Uii = U[i][j];  
		
			pi = computePressure(Ui); 
			pii = computePressure(Uii); 
			dp = fabs(pii - pi) / min(pii, pi);
			weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1)); 		

			Up = weight * Ui + (1 - weight) * Uii; 
			Um = (1 - weight) * Ui + weight * Uii; 

			Si = compute_inviscid_state(Up, nx, ny); 
			Sii = compute_inviscid_state(Um, nx, ny); 

			l = 0.5 * (Si.uprime + fabs(Si.uprime));
			lc = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - l); 
			lt = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a))); 

			V1_Plus = { lc * Si.k,
				(Si.u * lc + Si.a * nx * lt) * Si.k,
				(Si.v * lc + Si.a * ny * lt) * Si.k,
				(Si.h0 * lc + Si.a * Si.uprime * lt) * Si.k };

			V2_Plus = { lt / Si.a,
				Si.u * lt / Si.a + nx * lc,
				Si.v * lt / Si.a + ny * lc,
				Si.h0 * lt / Si.a + Si.uprime * lc };

			m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
			n = { -Si.uprime, nx, ny, 0 };

			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity();

			m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
			n = { -Sii.uprime, nx, ny, 0 };

			l = 0.5 * (Sii.uprime - fabs(Sii.uprime));
			lc = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - l);
			lt = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

			V1_Minus = { lc * Sii.k,
				(Sii.u * lc + Sii.a * nx * lt) * Sii.k,
				(Sii.v * lc + Sii.a * ny * lt) * Sii.k,
				(Sii.h0 * lc + Sii.a * Sii.uprime * lt) * Sii.k };

			V2_Minus = { lt / Sii.a,
				Sii.u * lt / Sii.a + nx * lc,
				Sii.v * lt / Sii.a + ny * lc,
				Sii.h0 * lt / Sii.a + Sii.uprime * lc };

			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity(); 

			j_plus_inviscid_Jacobians[i][j] = X; 
			j_minus_inviscid_Jacobians[i][j] = Y;
			j_Fluxes[i][j] = X * Ui + Y * Uii;   
		}
	}

} 

void Solver::compute_viscous_jacobians() {  

	static Vector Ui, Uii, Up, Um, V1_Plus, V2_Plus, V1_Minus, V2_Minus, n, m, dv, W;   
	Viscous_State Si, Sii; 
	double g = 5.72;
	double weight, dp, pi, pii, pe, nx, ny, lp, lcp, ltp, lm, lcm, ltm, rho, u, v, T, mu, lambda, k, dn, dx, dy;  
	pe = (gamma - 1);
	Matrix Mxx, Mxy, Myx, Myy, N, X, Y, Z, Elv, Erv, Ebv, Etv;   

	// Calculate Jacobians and Explicit fluxes for i-faces on left boundary 
	for (int j = 0; j < Ny; ++j) {

		Uii = U[0][j];

		Elv = viscous_boundary_2D_E(BoundaryType.left, Uii, grid.iNorms(0, j));   
		Ui = viscous_boundary_2D_U(BoundaryType.left, Uii, grid.iNorms(0, j)); 

		nx = grid.iNorms(0, j).x;
		ny = grid.iNorms(0, j).y;

		pi = computePressure(Ui);
		pii = computePressure(Uii);
		dp = fabs(pii - pi) / min(pii, pi);
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

		Up = weight * Ui + (1 - weight) * Uii;
		Um = (1 - weight) * Ui + weight * Uii;

		Si = compute_viscous_state(Up, nx, ny);
		Sii = compute_viscous_state(Um, nx, ny);

		W = 0.5 * (Ui + Uii);  
		rho = W[0]; 
		u = W[1] / W[0];  
		v = W[2] / W[0];
		T = computeTemperature(W); 
		mu = 1.458 * 1e-6 * T * sqrt(T) / (T + 110.3);
		lambda = -2 / 3 * mu;  
		k = cp * mu / Pr;

		Mxx = {
				{ {0, 0, 0, 0},
				  {0, lambda + 2 * mu, 0, 0},
				  {0, 0, mu, 0},
				  {0, u * (lambda + 2 * mu), v * mu, k}
				}
		};

		Mxy = {
			{ {0, 0, 0, 0},
			  {0, 0, lambda, 0},
			  {0, mu, 0, 0},
			  {0, v * mu, u * lambda, 0}
			}
		};

		Myx = {
			{ {0, 0, 0, 0},
			  {0, 0, mu, 0},
			  {0, lambda, 0, 0},
			  {0, v * lambda, u * mu, 0}
			}
		};

		Myy = {
				{ {0, 0, 0, 0},
				  {0, mu, 0, 0},
				  {0, 0, lambda + 2 * mu, 0},
				  {0, u * mu, v * (lambda + 2 * mu), k}
				}
		};

		N = {
			{ {1, 0, 0, 0},
			{-u / rho, 1 / rho, 0, 0},
			{-v / rho, 0, 1 / rho, 0},
			{ (u * u + v * v) / (2 * cv * rho) - T / rho, -u / (cv * rho), -v / (cv * rho), 1 / (cv * rho)}
			}
		};
	
		dx = grid.Center(1, j).x - grid.Center(0, j).x;
		dy = grid.Center(1, j).y - grid.Center(0, j).y; 
		cout << grid.Center(1, j).x << endl; 
		dn = sqrt(dx * dx + dy * dy); 

		dv = constoViscPrim(Uii) - constoViscPrim(Ui);
		Z = (Mxx * nx * nx + (Mxy + Myx) * nx * ny + Myy * ny * ny) / (-dn); 


		lp = 0.5 * (Si.uprime + fabs(Si.uprime)); 
		lcp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - lp);
		ltp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));

		V1_Plus = { lcp * Si.k,
			(Si.u * lcp + Si.a * nx * ltp) * Si.k,
			(Si.v * lcp + Si.a * ny * ltp) * Si.k,
			(Si.h0 * lcp + Si.a * Si.uprime * ltp) * Si.k };

		V2_Plus = { ltp / Si.a,
			Si.u * ltp / Si.a + nx * lcp,
			Si.v * ltp / Si.a + ny * lcp,
			Si.h0 * ltp / Si.a + Si.uprime * lcp };

		m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
		n = { -Si.uprime, nx, ny, 0 };

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity();

		m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
		n = { -Sii.uprime, nx, ny, 0 };

		lm = 0.5 * (Sii.uprime - fabs(Sii.uprime));
		lcm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - lm);
		ltm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

		V1_Minus = { lcm * Sii.k,
				(Sii.u * lcm + Sii.a * nx * ltm) * Sii.k,
				(Sii.v * lcm + Sii.a * ny * ltm) * Sii.k,
				(Sii.h0 * lcm + Sii.a * Sii.uprime * ltm) * Sii.k };

		V2_Minus = { ltm / Sii.a,
			Sii.u * ltm / Sii.a + nx * lcm,
			Sii.v * ltm / Sii.a + ny * lcm,
			Sii.h0 * ltm / Sii.a + Sii.uprime * lcm };


		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity();

		i_viscous_Jacobians[0][j] = Z * (identity() - Elv) * N; 

		//displayMatrix(Z);  
		i_plus_inviscid_Jacobians[0][j] = X;
		i_minus_inviscid_Jacobians[0][j] = Y;
		i_Fluxes[0][j] = X * Ui + Y * Uii + Z * dv;  
	}

	// Calculate Jacobians and Explicit fluxes for i-faces on right boundary
	for (int j = 0; j < Ny; ++j) {

		Ui = U[Nx - 1][j];
		Uii = viscous_boundary_2D_U(BoundaryType.right, Ui, grid.iNorms(Nx, j));
		Erv = viscous_boundary_2D_E(BoundaryType.right, Ui, grid.iNorms(Nx, j));

		nx = grid.iNorms(Nx, j).x;
		ny = grid.iNorms(Nx, j).y;

		pi = computePressure(Ui);
		pii = computePressure(Uii);
		dp = fabs(pii - pi) / min(pii, pi);
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

		Up = weight * Ui + (1 - weight) * Uii;
		Um = (1 - weight) * Ui + weight * Uii;

		Si = compute_viscous_state(Up, nx, ny);
		Sii = compute_viscous_state(Um, nx, ny);

		W = 0.5 * (Ui + Uii);
		rho = W[0];
		u = W[1] / W[0];
		v = W[2] / W[0];
		T = computeTemperature(W);
		mu = 1.458 * 1e-6 * T * sqrt(T) / (T + 110.3);
		lambda = -2 / 3 * mu;
		k = cp * mu / Pr;

		Mxx = {
				{ {0, 0, 0, 0},
				  {0, lambda + 2 * mu, 0, 0},
				  {0, 0, mu, 0},
				  {0, u * (lambda + 2 * mu), v * mu, k}
				}
		};

		Mxy = {
			{ {0, 0, 0, 0},
			  {0, 0, lambda, 0},
			  {0, mu, 0, 0},
			  {0, v * mu, u * lambda, 0}
			}
		};

		Myx = {
			{ {0, 0, 0, 0},
			  {0, 0, mu, 0},
			  {0, lambda, 0, 0},
			  {0, v * lambda, u * mu, 0}
			}
		};

		Myy = {
				{ {0, 0, 0, 0},
				  {0, mu, 0, 0},
				  {0, 0, lambda + 2 * mu, 0},
				  {0, u * mu, v * (lambda + 2 * mu), k}
				}
		};

		N = {
			{ {1, 0, 0, 0},
			{-u / rho, 1 / rho, 0, 0},
			{-v / rho, 0, 1 / rho, 0},
			{ (u * u + v * v) / (2 * cv * rho) - T / rho, -u / (cv * rho), -v / (cv * rho), 1 / (cv * rho)}
			}
		};

		dx = grid.Center(Nx - 1, j).x - grid.Center(Nx - 2, j).x; 
		dy = grid.Center(Nx - 1, j).y - grid.Center(Nx - 2, j).y;  

		dn = sqrt(dx * dx + dy * dy);
		dv = constoViscPrim(Uii) - constoViscPrim(Ui);

		Z = (Mxx * nx * nx + (Mxy + Myx) * nx * ny + Myy * ny * ny) / (-dn);

		lp = 0.5 * (Si.uprime + fabs(Si.uprime));
		lcp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - lp);
		ltp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));

		V1_Plus = { lcp * Si.k,
			(Si.u * lcp + Si.a * nx * ltp) * Si.k,
			(Si.v * lcp + Si.a * ny * ltp) * Si.k,
			(Si.h0 * lcp + Si.a * Si.uprime * ltp) * Si.k };

		V2_Plus = { ltp / Si.a,
			Si.u * ltp / Si.a + nx * lcp,
			Si.v * ltp / Si.a + ny * lcp,
			Si.h0 * ltp / Si.a + Si.uprime * lcp };

		m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
		n = { -Si.uprime, nx, ny, 0 };

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity();

		m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
		n = { -Sii.uprime, nx, ny, 0 };

		lm = 0.5 * (Sii.uprime - fabs(Sii.uprime));
		lcm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - lm);
		ltm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

		V1_Minus = { lcm * Sii.k,
				(Sii.u * lcm + Sii.a * nx * ltm) * Sii.k,
				(Sii.v * lcm + Sii.a * ny * ltm) * Sii.k,
				(Sii.h0 * lcm + Sii.a * Sii.uprime * ltm) * Sii.k };

		V2_Minus = { ltm / Sii.a,
			Sii.u * ltm / Sii.a + nx * lcm,
			Sii.v * ltm / Sii.a + ny * lcm,
			Sii.h0 * ltm / Sii.a + Sii.uprime * lcm };

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity();

		i_viscous_Jacobians[Nx][j] = Z * (Erv - identity()) * N;
		i_plus_inviscid_Jacobians[Nx][j] = X;
		i_minus_inviscid_Jacobians[Nx][j] = Y;
		i_Fluxes[Nx][j] = X * Ui + Y * Uii + Z * dv; 
	}

	// Calculate Jacobians and Explicit fluxes for j-faces on bottom boundary
	for (int i = 0; i < Nx; ++i) {

		Uii = U[i][0];
		Ui = viscous_boundary_2D_U(BoundaryType.bottom, Uii, grid.jNorms(i, 0));
		Ebv = viscous_boundary_2D_E(BoundaryType.bottom, Uii, grid.jNorms(i, 0)); 
		

		nx = grid.jNorms(i, 0).x;
		ny = grid.jNorms(i, 0).y;
	
		pi = computePressure(Ui);
		pii = computePressure(Uii);
		dp = fabs(pii - pi) / min(pii, pi);
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

		Up = weight * Ui + (1 - weight) * Uii;
		Um = (1 - weight) * Ui + weight * Uii;
		Si = compute_viscous_state(Up, nx, ny); 
		Sii = compute_viscous_state(Um, nx, ny); 

		W = 0.5 * (Ui + Uii);
		rho = W[0];
		u = W[1] / W[0];
		v = W[2] / W[0];
		T = computeTemperature(W);

		mu = 1.458 * 1e-6 * T * sqrt(T) / (T + 110.3);
		lambda = -2 / 3 * mu;
		k = cp * mu / Pr;

		Mxx = {
				{ {0, 0, 0, 0},
				  {0, lambda + 2 * mu, 0, 0},
				  {0, 0, mu, 0},
				  {0, u * (lambda + 2 * mu), v * mu, k}
				}
		};

		Mxy = {
			{ {0, 0, 0, 0},
			  {0, 0, lambda, 0},
			  {0, mu, 0, 0},
			  {0, v * mu, u * lambda, 0}
			}
		};

		Myx = {
			{ {0, 0, 0, 0},
			  {0, 0, mu, 0},
			  {0, lambda, 0, 0},
			  {0, v * lambda, u * mu, 0}
			}
		};

		Myy = {
				{ {0, 0, 0, 0},
				  {0, mu, 0, 0},
				  {0, 0, lambda + 2 * mu, 0},
				  {0, u * mu, v * (lambda + 2 * mu), k}
				}
		};

		N = {
			{ {1, 0, 0, 0},
			{-u / rho, 1 / rho, 0, 0},
			{-v / rho, 0, 1 / rho, 0},
			{ (u * u + v * v) / (2 * cv * rho) - T / rho, -u / (cv * rho), -v / (cv * rho), 1 / (cv * rho)}
			}
		};

		dx = grid.Center(i, 1).x - grid.Center(i, 0).x;   
		dy = grid.Center(i, 1).y - grid.Center(i, 0).y;    
		 
		dn = sqrt(dx * dx + dy * dy);
		dv = constoViscPrim(Uii) - constoViscPrim(Ui);

		Z = (Mxx * nx * nx + (Mxy + Myx) * nx * ny + Myy * ny * ny) / (-dn);

		lp = 0.5 * (Si.uprime + fabs(Si.uprime));
		lcp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - lp);
		ltp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));

		V1_Plus = { lcp * Si.k,
				(Si.u * lcp + Si.a * nx * ltp) * Si.k,
				(Si.v * lcp + Si.a * ny * ltp) * Si.k,
				(Si.h0 * lcp + Si.a * Si.uprime * ltp) * Si.k };

		V2_Plus = { ltp / Si.a,
			Si.u * ltp / Si.a + nx * lcp,
			Si.v * ltp / Si.a + ny * lcp,
			Si.h0 * ltp / Si.a + Si.uprime * lcp };

		m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
		n = { -Si.uprime, nx, ny, 0 };

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity();

		m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
		n = { -Sii.uprime, nx, ny, 0 };

		lm = 0.5 * (Sii.uprime - fabs(Sii.uprime));
		lcm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - lm);
		ltm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

		V1_Minus = { lcm * Sii.k,
				(Sii.u * lcm + Sii.a * nx * ltm) * Sii.k,
				(Sii.v * lcm + Sii.a * ny * ltm) * Sii.k,
				(Sii.h0 * lcm + Sii.a * Sii.uprime * ltm) * Sii.k };

		V2_Minus = { ltm / Sii.a,
			Sii.u * ltm / Sii.a + nx * lcm,
			Sii.v * ltm / Sii.a + ny * lcm,
			Sii.h0 * ltm / Sii.a + Sii.uprime * lcm };

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity();

		j_viscous_Jacobians[i][0] = Z * (identity() - Ebv) * N;  
		j_plus_inviscid_Jacobians[i][0] = X;
		j_minus_inviscid_Jacobians[i][0] = Y;
		j_Fluxes[i][0] = X * Ui + Y * Uii + Z * dv; 

	}

	// Calculate Jacobians and Explicit fluxes for j-faces on top boundary
	for (int i = 0; i < Nx; ++i) {

		Ui = U[i][Ny - 1];
		Uii = viscous_boundary_2D_U(BoundaryType.top, Ui, grid.jNorms(i, Ny));
		Etv = viscous_boundary_2D_E(BoundaryType.top, Ui, grid.jNorms(i, Ny));

		nx = grid.jNorms(i, Ny).x;
		ny = grid.jNorms(i, Ny).y;

		pi = computePressure(Ui);
		pii = computePressure(Uii);
		dp = fabs(pii - pi) / min(pii, pi);
		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

		Up = weight * Ui + (1 - weight) * Uii;
		Um = (1 - weight) * Ui + weight * Uii;

		Si = compute_viscous_state(Up, nx, ny);
		Sii = compute_viscous_state(Um, nx, ny);

		W = 0.5 * (Ui + Uii);
		rho = W[0];
		u = W[1] / W[0];
		v = W[2] / W[0];
		T = computeTemperature(W);
		mu = 1.458 * 1e-6 * T * sqrt(T) / (T + 110.3);
		lambda = -2 / 3 * mu;
		k = cp * mu / Pr;

		Mxx = {
				{ {0, 0, 0, 0},
				  {0, lambda + 2 * mu, 0, 0},
				  {0, 0, mu, 0},
				  {0, u * (lambda + 2 * mu), v * mu, k}
				}
		};

		Mxy = {
			{ {0, 0, 0, 0},
			  {0, 0, lambda, 0},
			  {0, mu, 0, 0},
			  {0, v * mu, u * lambda, 0}
			}
		};

		Myx = {
			{ {0, 0, 0, 0},
			  {0, 0, mu, 0},
			  {0, lambda, 0, 0},
			  {0, v * lambda, u * mu, 0}
			}
		};

		Myy = {
				{ {0, 0, 0, 0},
				  {0, mu, 0, 0},
				  {0, 0, lambda + 2 * mu, 0},
				  {0, u * mu, v * (lambda + 2 * mu), k}
				}
		};

		N = {
			{ {1, 0, 0, 0},
			{-u / rho, 1 / rho, 0, 0},
			{-v / rho, 0, 1 / rho, 0},
			{ (u * u + v * v) / (2 * cv * rho) - T / rho, -u / (cv * rho), -v / (cv * rho), 1 / (cv * rho)}
			}
		};

		dx = grid.Center(i, Ny - 1).x - grid.Center(i, Ny - 2).x;
		dy = grid.Center(i, Ny - 1).y - grid.Center(i, Ny - 2).y;   

		dn = sqrt(dx * dx + dy * dy);
		dv = constoViscPrim(Uii) - constoViscPrim(Ui);

		Z = (Mxx * nx * nx + (Mxy + Myx) * nx * ny + Myy * ny * ny) / (-dn); 

		lp = 0.5 * (Si.uprime + fabs(Si.uprime));
		lcp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - lp);
		ltp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));

		V1_Plus = { lcp * Si.k,
				(Si.u * lcp + Si.a * nx * ltp) * Si.k,
				(Si.v * lcp + Si.a * ny * ltp) * Si.k,
				(Si.h0 * lcp + Si.a * Si.uprime * ltp) * Si.k };

		V2_Plus = { ltp / Si.a,
			Si.u * ltp / Si.a + nx * lcp,
			Si.v * ltp / Si.a + ny * lcp,
			Si.h0 * ltp / Si.a + Si.uprime * lcp };

		m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
		n = { -Si.uprime, nx, ny, 0 };

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity();

		m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
		n = { -Sii.uprime, nx, ny, 0 };

		lm = 0.5 * (Sii.uprime - fabs(Sii.uprime));
		lcm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - lm);
		ltm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

		V1_Minus = { lcm * Sii.k,
				(Sii.u * lcm + Sii.a * nx * ltm) * Sii.k,
				(Sii.v * lcm + Sii.a * ny * ltm) * Sii.k,
				(Sii.h0 * lcm + Sii.a * Sii.uprime * ltm) * Sii.k };

		V2_Minus = { ltm / Sii.a,
			Sii.u * ltm / Sii.a + nx * lcm,
			Sii.v * ltm / Sii.a + ny * lcm,
			Sii.h0 * ltm / Sii.a + Sii.uprime * lcm };

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity();

		j_viscous_Jacobians[i][Ny] = Z * (Etv - identity()) * N; 
		j_plus_inviscid_Jacobians[i][Ny] = X;
		j_minus_inviscid_Jacobians[i][Ny] = Y;
		j_Fluxes[i][Ny] = X * Ui + Y * Uii + Z * dv; 

	}

	// Inner i-faces  
	for (int i = 1; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {

			nx = grid.iNorms(i, j).x;
			ny = grid.iNorms(i, j).y;

			Ui = U[i - 1][j];
			Uii = U[i][j];

			pi = computePressure(Ui);
			pii = computePressure(Uii);
			dp = fabs(pii - pi) / min(pii, pi);
			weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));


			Up = weight * Ui + (1 - weight) * Uii;
			Um = (1 - weight) * Ui + weight * Uii;

			Si = compute_viscous_state(Up, nx, ny);
			Sii = compute_viscous_state(Um, nx, ny);

			W = 0.5 * (Ui + Uii);
			rho = W[0];
			u = W[1] / W[0];
			v = W[2] / W[0];
			T = computeTemperature(W);
			mu = 1.458 * 1e-6 * T * sqrt(T) / (T + 110.3);
			lambda = -2 / 3 * mu;
			k = cp * mu / Pr;

			Mxx = {
				{ {0, 0, 0, 0},
				  {0, lambda + 2*mu, 0, 0},
				  {0, 0, mu, 0},
				  {0, u * (lambda + 2 * mu), v * mu, k}
				}
			};

			Mxy = {
				{ {0, 0, 0, 0},
				  {0, 0, lambda, 0},
				  {0, mu, 0, 0},
				  {0, v * mu, u * lambda, 0}  
				}
			};

			Myx = {
				{ {0, 0, 0, 0},
				  {0, 0, mu, 0},
				  {0, lambda, 0, 0},
				  {0, v * lambda, u * mu, 0}
				}
			};

			Myy = {
				{ {0, 0, 0, 0},
				  {0, mu, 0, 0},
				  {0, 0, lambda + 2 * mu, 0},
				  {0, u * mu, v * (lambda + 2 * mu), k}
				}
			};

			N = {
				{ {1, 0, 0, 0},
				{-u / rho, 1 / rho, 0, 0},
				{-v / rho, 0, 1 / rho, 0},
				{ (u * u + v * v) / (2 * cv * rho) - T / rho, -u / (cv * rho), -v / (cv * rho), 1 / (cv * rho)}
				}
			};

			dx = grid.Center(i, j).x - grid.Center(i - 1, j).x;  
			dy = grid.Center(i, j).y - grid.Center(i - 1, j).y;  

			dn = sqrt(dx * dx + dy * dy);
			dv = constoViscPrim(Uii) - constoViscPrim(Ui);

			Z = (Mxx * nx * nx + (Mxy + Myx) * nx * ny + Myy * ny * ny) / (-dn);

			lp = 0.5 * (Si.uprime + fabs(Si.uprime));
			lcp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - lp);
			ltp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));

			V1_Plus = { lcp * Si.k,
				(Si.u * lcp + Si.a * nx * ltp) * Si.k,
				(Si.v * lcp + Si.a * ny * ltp) * Si.k,
				(Si.h0 * lcp + Si.a * Si.uprime * ltp) * Si.k };

			V2_Plus = { ltp / Si.a,
				Si.u * ltp / Si.a + nx * lcp,
				Si.v * ltp / Si.a + ny * lcp,
				Si.h0 * ltp / Si.a + Si.uprime * lcp };

			m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
			n = { -Si.uprime, nx, ny, 0 };

			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity();


			m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
			n = { -Sii.uprime, nx, ny, 0 };

			lm = 0.5 * (Sii.uprime - fabs(Sii.uprime));
			lcm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - lm);
			ltm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

			V1_Minus = { lcm * Sii.k,
				(Sii.u * lcm + Sii.a * nx * ltm) * Sii.k,
				(Sii.v * lcm + Sii.a * ny * ltm) * Sii.k,
				(Sii.h0 * lcm + Sii.a * Sii.uprime * ltm) * Sii.k };

			V2_Minus = { ltm / Sii.a,
				Sii.u * ltm / Sii.a + nx * lcm,
				Sii.v * ltm / Sii.a + ny * lcm,
				Sii.h0 * ltm / Sii.a + Sii.uprime * lcm };

			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity();

			i_viscous_Jacobians[i][j] = Z * N;  
			i_plus_inviscid_Jacobians[i][j] = X;
			i_minus_inviscid_Jacobians[i][j] = Y;
			i_Fluxes[i][j] = X * Ui + Y * Uii + Z * dv;  
		}
	}

	// Inner j-faces 
	for (int i = 0; i < Nx; ++i) {
		for (int j = 1; j < Ny; ++j) {

			nx = grid.jNorms(i, j).x;
			ny = grid.jNorms(i, j).y;

			Ui = U[i][j - 1];
			Uii = U[i][j];

			pi = computePressure(Ui);
			pii = computePressure(Uii);
			dp = fabs(pii - pi) / min(pii, pi);
			weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));

			Up = weight * Ui + (1 - weight) * Uii;
			Um = (1 - weight) * Ui + weight * Uii;

			Si = compute_viscous_state(Up, nx, ny);
			Sii = compute_viscous_state(Um, nx, ny);

			W = 0.5 * (Ui + Uii);
			rho = W[0];
			u = W[1] / W[0];
			v = W[2] / W[0];
			T = computeTemperature(W);
			mu = 1.458 * 1e-6 * T * sqrt(T) / (T + 110.3);
			lambda = -2 / 3 * mu;
			k = cp * mu / Pr;

			Mxx = {
				{ {0, 0, 0, 0},
				  {0, lambda + 2 * mu, 0, 0},
				  {0, 0, mu, 0},
				  {0, u * (lambda + 2 * mu), v * mu, k}
				}
			};

			Mxy = {
				{ {0, 0, 0, 0},
				  {0, 0, lambda, 0},
				  {0, mu, 0, 0},
				  {0, v * mu, u * lambda, 0}
				}
			};

			Myx = {
				{ {0, 0, 0, 0},
				  {0, 0, mu, 0},
				  {0, lambda, 0, 0},
				  {0, v * lambda, u * mu, 0}
				}
			};

			Myy = {
				{ {0, 0, 0, 0},
				  {0, mu, 0, 0},
				  {0, 0, lambda + 2 * mu, 0},
				  {0, u * mu, v * (lambda + 2 * mu), k}
				}
			};


			N = {
				{ {1, 0, 0, 0},
				{-u / rho, 1 / rho, 0, 0},
				{-v / rho, 0, 1 / rho, 0},
				{ (u * u + v * v) / (2 * cv * rho) - T / rho, -u / (cv * rho), -v / (cv * rho), 1 / (cv * rho)}
				}
			};

			dx = grid.Center(i, j).x - grid.Center(i, j - 1).x;
			dy = grid.Center(i, j).y - grid.Center(i, j - 1).y; 

			dn = sqrt(dx * dx + dy * dy);
			dv = constoViscPrim(Uii) - constoViscPrim(Ui);

			Z = (Mxx * nx * nx + (Mxy + Myx) * nx * ny + Myy * ny * ny) / (-dn);

			lp = 0.5 * (Si.uprime + fabs(Si.uprime));
			lcp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - lp);
			ltp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));

			V1_Plus = { lcp * Si.k,
				(Si.u * lcp + Si.a * nx * ltp) * Si.k,
				(Si.v * lcp + Si.a * ny * ltp) * Si.k,
				(Si.h0 * lcp + Si.a * Si.uprime * ltp) * Si.k };

			V2_Plus = { ltp / Si.a,
				Si.u * ltp / Si.a + nx * lcp,
				Si.v * ltp / Si.a + ny * lcp,
				Si.h0 * ltp / Si.a + Si.uprime * lcp };

			m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
			n = { -Si.uprime, nx, ny, 0 };

			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity();

			m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
			n = { -Sii.uprime, nx, ny, 0 };

			lm = 0.5 * (Sii.uprime - fabs(Sii.uprime));
			lcm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - lm);
			ltm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

			V1_Minus = { lcm * Sii.k,
				(Sii.u * lcm + Sii.a * nx * ltm) * Sii.k,
				(Sii.v * lcm + Sii.a * ny * ltm) * Sii.k,
				(Sii.h0 * lcm + Sii.a * Sii.uprime * ltm) * Sii.k };

			V2_Minus = { ltm / Sii.a,
				Sii.u * ltm / Sii.a + nx * lcm,
				Sii.v * ltm / Sii.a + ny * lcm,
				Sii.h0 * ltm / Sii.a + Sii.uprime * lcm };

			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity();

			j_viscous_Jacobians[i][j] = Z * N; 
			j_plus_inviscid_Jacobians[i][j] = X;
			j_minus_inviscid_Jacobians[i][j] = Y;
			j_Fluxes[i][j] = X * Ui + Y * Uii + Z * dv; 
		}
	}

}

void Solver::solve_left_line_inviscid() {

	int i = 0; 
	Matrix Et, Eb, El; 
	static Matrix A; // Relates to U(j+1) 
	static Matrix B; // Relates to U(j) 
	static Matrix C; // Relates to U(j-1) 
	static Vector F; // Right hand side  

	// Intermediate matrices
	static Matrix alpha;
	static array<array<double, 4>, Ny> v;
	static array<array<array<double,4>, 4>, Ny> g; 

	// Grab boundary values
	Eb = inviscid_boundary_2D_E(BoundaryType.bottom, U[i][0], grid.jNorms(i, 0));      
	Et = inviscid_boundary_2D_E(BoundaryType.top, U[i][Ny - 1], grid.jNorms(i, Ny));     
	El = inviscid_boundary_2D_E(BoundaryType.left, U[i][Ny - 1], grid.iNorms(i, Ny - 1)); 

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity()
		- (i_minus_inviscid_Jacobians[i][Ny - 1] + El * i_plus_inviscid_Jacobians[i][Ny - 1]) * grid.iArea(i, Ny - 1)	// Left 
		+ (i_plus_inviscid_Jacobians[i + 1][Ny - 1]) * grid.iArea(i + 1, Ny - 1)										// Right
		- (j_minus_inviscid_Jacobians[i][Ny - 1]) * grid.jArea(i, Ny - 1)												// Bottom
		+ (j_plus_inviscid_Jacobians[i][Ny] + Et * j_minus_inviscid_Jacobians[i][Ny]) * grid.jArea(i, Ny);				// Top 

	C = j_plus_inviscid_Jacobians[i][Ny - 1] * (-grid.jArea(i, Ny - 1));

	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1))
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny))
		+ i_minus_inviscid_Jacobians[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1)) * dU_old[i + 1][Ny - 1];

	alpha = A;
	v[Ny - 1] = F / alpha;
	g[Ny - 1] = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		El = inviscid_boundary_2D_E(BoundaryType.left, U[i][j], grid.iNorms(i, j));  

		B = j_minus_inviscid_Jacobians[i][j + 1] * (grid.jArea(i, j + 1));

		A = grid.Volume(i, j) / dt * identity()
			- (i_minus_inviscid_Jacobians[i][j] + El * i_plus_inviscid_Jacobians[i][j]) * grid.iArea(i, j)	// Left
			+ i_plus_inviscid_Jacobians[i + 1][j] * grid.iArea(i + 1, j)									// Right
			- j_minus_inviscid_Jacobians[i][j] * grid.jArea(i, j)											// Bottom
			+ j_plus_inviscid_Jacobians[i][j + 1] * grid.jArea(i, j + 1);									// Top

		C = j_plus_inviscid_Jacobians[i][j] * (-grid.jArea(i, j));

		F = i_Fluxes[i][j] * (grid.iArea(i, j))
			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j))
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1))
			+ i_minus_inviscid_Jacobians[i + 1][j] * (-grid.iArea(i + 1, j)) * dU_old[i + 1][j];



		alpha = A - B * g[j + 1];
		g[j] = C / alpha;
		v[j] = (F - B * v[j + 1]) / alpha;
	}

	//  Bottom boundary

	El = inviscid_boundary_2D_E(BoundaryType.left, U[i][0], grid.iNorms(i, 0)); 

	B = j_minus_inviscid_Jacobians[i][1] * (grid.jArea(i, 1));

	A = grid.Volume(i, 0) / dt * identity()
		- (i_minus_inviscid_Jacobians[i][0] + El * i_plus_inviscid_Jacobians[i][0]) * grid.iArea(i, 0)		// Left
		+ i_plus_inviscid_Jacobians[i + 1][0] * grid.iArea(i + 1, 0)										// Right
		- (Eb * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0]) * grid.jArea(i, 0)		// Bottom
		+ j_plus_inviscid_Jacobians[i][1] * grid.jArea(i, 1);												// Top

	F = i_Fluxes[i][0] * (grid.iArea(i, 0))
		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1))
		+ i_minus_inviscid_Jacobians[i + 1][0] * (-grid.iArea(i + 1, 0)) * dU_old[i + 1][0];


	alpha = A - B * g[1];
	v[0] = (F - B * v[1]) / alpha;

	// Calculate dU
	dU_new[i][0] = v[0];
	for (int j = 1; j < Ny; ++j) {
		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];
	}

}

void Solver::solve_middle_line_inviscid(const int i) {

	//auto start = TIME;
	Matrix Et, Eb; 
	static Matrix A; // Relates to U(j+1) 
	static Matrix B; // Relates to U(j) 
	static Matrix C; // Relates to U(j-1) 
	static Vector F; // Right hand side  

	// Intermediate matrices
	static Matrix alpha; 
	static array<array<double, 4>, Ny> v;
	static array<array<array<double, 4>, 4>, Ny> g; 

	// Grab boundary values
	Eb = inviscid_boundary_2D_E(BoundaryType.bottom, U[i][0], grid.jNorms(i, 0));  
	Et = inviscid_boundary_2D_E(BoundaryType.top, U[i][Ny - 1], grid.jNorms(i, Ny));    

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity() 
		- i_minus_inviscid_Jacobians[i][Ny - 1] * grid.iArea(i, Ny - 1) 
		+ i_plus_inviscid_Jacobians[i + 1][Ny - 1] * grid.iArea(i + 1, Ny - 1)
		- j_minus_inviscid_Jacobians[i][Ny - 1] * grid.jArea(i, Ny - 1) 
		+ (j_plus_inviscid_Jacobians[i][Ny] + Et * j_minus_inviscid_Jacobians[i][Ny]) * grid.jArea(i, Ny); 

	C = j_plus_inviscid_Jacobians[i][Ny - 1] * (-grid.jArea(i, Ny - 1)); 

	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1)) 
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny))
		+ i_plus_inviscid_Jacobians[i][Ny - 1] * (grid.iArea(i, Ny - 1)) * dU_old[i - 1][Ny - 1]
		+ i_minus_inviscid_Jacobians[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1)) * dU_old[i + 1][Ny - 1]; 
	
	alpha = A;
	v[Ny - 1] = F / alpha;  
	g[Ny - 1] = C / alpha; 

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		B = j_minus_inviscid_Jacobians[i][j + 1] * (grid.jArea(i, j + 1)); 

		A = grid.Volume(i, j) / dt * identity() 
			- i_minus_inviscid_Jacobians[i][j] * grid.iArea(i, j) 
			+ i_plus_inviscid_Jacobians[i + 1][j] * grid.iArea(i + 1, j)  
			- j_minus_inviscid_Jacobians[i][j] * grid.jArea(i, j)  
			+ j_plus_inviscid_Jacobians[i][j + 1] * grid.jArea(i, j + 1);   

		C = j_plus_inviscid_Jacobians[i][j] * (-grid.jArea(i, j));  

		F = i_Fluxes[i][j] * (grid.iArea(i, j))  
			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j)) 
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1)) 
			+ i_plus_inviscid_Jacobians[i][j] * (grid.iArea(i, j)) * dU_old[i - 1][j] 
			+ i_minus_inviscid_Jacobians[i + 1][j] * (-grid.iArea(i + 1, j)) * dU_old[i + 1][j]; 



		alpha = A - B * g[j + 1]; 
		g[j] = C / alpha; 
		v[j] = (F - B * v[j + 1]) / alpha; 
	}

	//  Bottom boundary
	B = j_minus_inviscid_Jacobians[i][1] * (grid.jArea(i, 1));

	A = grid.Volume(i, 0) / dt * identity()
		- i_minus_inviscid_Jacobians[i][0] * grid.iArea(i, 0)
		+ i_plus_inviscid_Jacobians[i + 1][0] * grid.iArea(i + 1, 0)
		- (Eb * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0]) * grid.jArea(i, 0) 
		+ j_plus_inviscid_Jacobians[i][1] * grid.jArea(i, 1); 

	F = i_Fluxes[i][0] * (grid.iArea(i, 0))
		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1))
		+ i_plus_inviscid_Jacobians[i][0] * (grid.iArea(i, 0)) * dU_old[i - 1][0]
		+ i_minus_inviscid_Jacobians[i + 1][0] * (-grid.iArea(i + 1, 0)) * dU_old[i + 1][0];


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

void Solver::solve_right_line_inviscid() {

	int i = Nx - 1; 
	Matrix Et, Eb, Er;  

	static Matrix A; // Relates to U(j+1) 
	static Matrix B; // Relates to U(j) 
	static Matrix C; // Relates to U(j-1) 
	static Vector F; // Right hand side  

	// Intermediate matrices
	static Matrix alpha;
	static array<array<double, 4>, Ny> v; 
	static array<array<array<double, 4>, 4>, Ny> g; 

	// Grab boundary values
	Eb = inviscid_boundary_2D_E(BoundaryType.bottom, U[i][0], grid.jNorms(i, 0));
	Et = inviscid_boundary_2D_E(BoundaryType.top, U[i][Ny - 1], grid.jNorms(i, Ny));
	Er = inviscid_boundary_2D_E(BoundaryType.right, U[i][Ny - 1], grid.iNorms(i, Ny - 1)); 

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity()
		- i_minus_inviscid_Jacobians[i][Ny - 1] * grid.iArea(i, Ny - 1)																// Left
		+ (i_plus_inviscid_Jacobians[i + 1][Ny - 1] + Er * i_minus_inviscid_Jacobians[i + 1][Ny - 1]) * grid.iArea(i + 1, Ny - 1)	// Right
		- j_minus_inviscid_Jacobians[i][Ny - 1] * grid.jArea(i, Ny - 1)																// Bottom
		+ (j_plus_inviscid_Jacobians[i][Ny] + Et * j_minus_inviscid_Jacobians[i][Ny]) * grid.jArea(i, Ny);							// Top 

	C = j_plus_inviscid_Jacobians[i][Ny - 1] * (-grid.jArea(i, Ny - 1));

	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1))
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny))
		+ i_plus_inviscid_Jacobians[i][Ny - 1] * (grid.iArea(i, Ny - 1)) * dU_old[i - 1][Ny - 1];

	alpha = A;
	v[Ny - 1] = F / alpha;
	g[Ny - 1] = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		Er = inviscid_boundary_2D_E(BoundaryType.right, U[i][j], grid.iNorms(i, j));  

		B = j_minus_inviscid_Jacobians[i][j + 1] * (grid.jArea(i, j + 1));

		A = grid.Volume(i, j) / dt * identity()
			- i_minus_inviscid_Jacobians[i][j] * grid.iArea(i, j)														// Left
			+ (i_plus_inviscid_Jacobians[i + 1][j] + Er * i_minus_inviscid_Jacobians[i + 1][j])* grid.iArea(i + 1, j)	// Right
			- j_minus_inviscid_Jacobians[i][j] * grid.jArea(i, j)														// Bottom
			+ j_plus_inviscid_Jacobians[i][j + 1] * grid.jArea(i, j + 1);												// Top

		C = j_plus_inviscid_Jacobians[i][j] * (-grid.jArea(i, j));

		F = i_Fluxes[i][j] * (grid.iArea(i, j))
			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j))
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1))
			+ i_plus_inviscid_Jacobians[i][j] * (grid.iArea(i, j)) * dU_old[i - 1][j];



		alpha = A - B * g[j + 1];
		g[j] = C / alpha;
		v[j] = (F - B * v[j + 1]) / alpha;
	}
	//  Bottom boundary
	Er = inviscid_boundary_2D_E(BoundaryType.right, U[i][0], grid.iNorms(i, 0)); 

	B = j_minus_inviscid_Jacobians[i][1] * (grid.jArea(i, 1));

	A = grid.Volume(i, 0) / dt * identity()
		- i_minus_inviscid_Jacobians[i][0] * grid.iArea(i, 0)															// Left
		+ (i_plus_inviscid_Jacobians[i + 1][0] + Er * i_minus_inviscid_Jacobians[i + 1][0]) * grid.iArea(i + 1, 0)		// Right
		- (Eb * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0]) * grid.jArea(i, 0)					// Bottom
		+ j_plus_inviscid_Jacobians[i][1] * grid.jArea(i, 1);															// Top 

	F = i_Fluxes[i][0] * (grid.iArea(i, 0))
		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1))
		+ i_plus_inviscid_Jacobians[i][0] * (grid.iArea(i, 0)) * dU_old[i - 1][0];


	alpha = A - B * g[1];
	v[0] = (F - B * v[1]) / alpha;

	// Calculate dU
	dU_new[i][0] = v[0];

	for (int j = 1; j < Ny; ++j) {
		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];
	}

	

}

void Solver::solve_left_line_viscous() {

	int i = 0;
	Matrix Eti, Ebi, Eli; 
	static Matrix A; // Relates to U(j+1) 
	static Matrix B; // Relates to U(j) 
	static Matrix C; // Relates to U(j-1) 
	static Vector F; // Right hand side  

	// Intermediate matrices
	static Matrix alpha;
	static array<array<double, 4>, Ny> v;
	static array<array<array<double, 4>, 4>, Ny> g;

	// Grab boundary values
	Ebi = inviscid_boundary_2D_E(BoundaryType.bottom, U[i][0], grid.jNorms(i, 0));
	Eti = inviscid_boundary_2D_E(BoundaryType.top, U[i][Ny - 1], grid.jNorms(i, Ny));

	Eli = inviscid_boundary_2D_E(BoundaryType.left, U[i][Ny - 1], grid.iNorms(i, Ny - 1));

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity()
		- (i_minus_inviscid_Jacobians[i][Ny - 1] + Eli * i_plus_inviscid_Jacobians[i][Ny - 1] + i_viscous_Jacobians[i][Ny - 1]) * grid.iArea(i, Ny - 1)		// Left     
		+ (i_plus_inviscid_Jacobians[i + 1][Ny - 1] - i_viscous_Jacobians[i + 1][Ny - 1]) * grid.iArea(i + 1, Ny - 1)										// Right
		- (j_minus_inviscid_Jacobians[i][Ny - 1] + j_viscous_Jacobians[i][Ny - 1]) * grid.jArea(i, Ny - 1)		 											// Bottom
		+ (j_plus_inviscid_Jacobians[i][Ny] + Eti * j_minus_inviscid_Jacobians[i][Ny] + j_viscous_Jacobians[i][Ny]) * grid.jArea(i, Ny);					// Top  

	C = (j_plus_inviscid_Jacobians[i][Ny - 1] - j_viscous_Jacobians[i][Ny - 1]) * (-grid.jArea(i, Ny - 1)); 

	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1))
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny))
		+ (i_minus_inviscid_Jacobians[i + 1][Ny - 1] + i_viscous_Jacobians[i + 1][Ny - 1]) * (-grid.iArea(i + 1, Ny - 1)) * dU_old[i + 1][Ny - 1];   

	alpha = A;
	v[Ny - 1] = F / alpha;
	g[Ny - 1] = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		Eli = inviscid_boundary_2D_E(BoundaryType.left, U[i][j], grid.iNorms(i, j)); 

		B = (j_minus_inviscid_Jacobians[i][j + 1] + j_viscous_Jacobians[i][j + 1]) * (grid.jArea(i, j + 1));

		A = grid.Volume(i, j) / dt * identity()
			- (i_minus_inviscid_Jacobians[i][j] + Eli * i_plus_inviscid_Jacobians[i][j] + i_viscous_Jacobians[i][j]) * grid.iArea(i, j)			// Left  
			+ (i_plus_inviscid_Jacobians[i + 1][j] - i_viscous_Jacobians[i + 1][j] ) * grid.iArea(i + 1, j)										// Right
			- (j_minus_inviscid_Jacobians[i][j] + j_viscous_Jacobians[i][j]) * grid.jArea(i, j)													// Bottom
			+ (j_plus_inviscid_Jacobians[i][j + 1] - j_viscous_Jacobians[i][j + 1]) * grid.jArea(i, j + 1);										// Top

		C = (j_plus_inviscid_Jacobians[i][j] - j_viscous_Jacobians[i][j]) * (-grid.jArea(i, j));  

		F = i_Fluxes[i][j] * (grid.iArea(i, j))
			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j))
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1))
			+ (i_minus_inviscid_Jacobians[i + 1][j] + i_viscous_Jacobians[i + 1][j]) * (-grid.iArea(i + 1, j)) * dU_old[i + 1][j]; 



		alpha = A - B * g[j + 1];
		g[j] = C / alpha;
		v[j] = (F - B * v[j + 1]) / alpha;
	}

	//  Bottom boundary

	Eli = inviscid_boundary_2D_E(BoundaryType.left, U[i][0], grid.iNorms(i, 0));

	B = (j_minus_inviscid_Jacobians[i][1] + j_viscous_Jacobians[i][1]) * (grid.jArea(i, 1)); 

	A = grid.Volume(i, 0) / dt * identity()
		- (i_minus_inviscid_Jacobians[i][0] + Eli * i_plus_inviscid_Jacobians[i][0] + i_viscous_Jacobians[i][0]) * grid.iArea(i, 0)		// Left
		+ (i_plus_inviscid_Jacobians[i + 1][0] - i_viscous_Jacobians[i + 1][0]) * grid.iArea(i + 1, 0)									// Right 
		- (Ebi * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0] + j_viscous_Jacobians[i][0]) * grid.jArea(i, 0)		// Bottom
		+ (j_plus_inviscid_Jacobians[i][1] - j_viscous_Jacobians[i][1]) * grid.jArea(i, 1);												// Top

	F = i_Fluxes[i][0] * (grid.iArea(i, 0))
		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1))
		+ (i_minus_inviscid_Jacobians[i + 1][0] + i_viscous_Jacobians[i + 1][0]) * (-grid.iArea(i + 1, 0)) * dU_old[i + 1][0]; 


	alpha = A - B * g[1];
	v[0] = (F - B * v[1]) / alpha;

	// Calculate dU
	dU_new[i][0] = v[0];
	for (int j = 1; j < Ny; ++j) {
		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];
	}

}

void Solver::solve_middle_line_viscous(const int i) {

	Matrix Eti, Ebi; 

	static Matrix A; // Relates to U(j+1) 
	static Matrix B; // Relates to U(j) 
	static Matrix C; // Relates to U(j-1) 
	static Vector F; // Right hand side  

	// Intermediate matrices
	static Matrix alpha;
	static array<array<double, 4>, Ny> v;
	static array<array<array<double, 4>, 4>, Ny> g;

	// Grab boundary values
	Ebi = inviscid_boundary_2D_E(BoundaryType.bottom, U[i][0], grid.jNorms(i, 0));
	Eti = inviscid_boundary_2D_E(BoundaryType.top, U[i][Ny - 1], grid.jNorms(i, Ny));
	

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity()
		- (i_minus_inviscid_Jacobians[i][Ny - 1] + i_viscous_Jacobians[i][Ny - 1]) * grid.iArea(i, Ny - 1)															// Left
		+ (i_plus_inviscid_Jacobians[i + 1][Ny - 1] - i_viscous_Jacobians[i + 1][Ny - 1]) * grid.iArea(i + 1, Ny - 1)												// Right
		- (j_minus_inviscid_Jacobians[i][Ny - 1] + j_viscous_Jacobians[i][Ny - 1]) * grid.jArea(i, Ny - 1)															// Bottom
		+ (j_plus_inviscid_Jacobians[i][Ny] + Eti * j_minus_inviscid_Jacobians[i][Ny] + j_viscous_Jacobians[i][Ny]) * grid.jArea(i, Ny);		// Right 

	C = (j_plus_inviscid_Jacobians[i][Ny - 1] - j_viscous_Jacobians[i][Ny - 1]) * (-grid.jArea(i, Ny - 1)); 

	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1))
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny))
		+ (i_plus_inviscid_Jacobians[i][Ny - 1] - i_viscous_Jacobians[i][Ny - 1]) * (grid.iArea(i, Ny - 1)) * dU_old[i - 1][Ny - 1]
		+ (i_minus_inviscid_Jacobians[i + 1][Ny - 1] + i_viscous_Jacobians[i + 1][Ny - 1]) * (-grid.iArea(i + 1, Ny - 1)) * dU_old[i + 1][Ny - 1]; 

	alpha = A;
	v[Ny - 1] = F / alpha;
	g[Ny - 1] = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		B = (j_minus_inviscid_Jacobians[i][j + 1] + j_viscous_Jacobians[i][j + 1]) * (grid.jArea(i, j + 1));

		A = grid.Volume(i, j) / dt * identity()
			- (i_minus_inviscid_Jacobians[i][j] + i_viscous_Jacobians[i][j]) * grid.iArea(i, j)
			+ (i_plus_inviscid_Jacobians[i + 1][j] - i_viscous_Jacobians[i + 1][j]) * grid.iArea(i + 1, j)
			- (j_minus_inviscid_Jacobians[i][j] + j_viscous_Jacobians[i][j]) * grid.jArea(i, j)
			+ (j_plus_inviscid_Jacobians[i][j + 1] - j_viscous_Jacobians[i][j + 1]) * grid.jArea(i, j + 1); 

		C = (j_plus_inviscid_Jacobians[i][j] - j_viscous_Jacobians[i][j]) * (-grid.jArea(i, j)); 

		F = i_Fluxes[i][j] * (grid.iArea(i, j))
			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j))
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1))
			+ (i_plus_inviscid_Jacobians[i][j] - i_viscous_Jacobians[i][j])* (grid.iArea(i, j)) * dU_old[i - 1][j]
			+ (i_minus_inviscid_Jacobians[i + 1][j] + i_viscous_Jacobians[i + 1][j]) * (-grid.iArea(i + 1, j)) * dU_old[i + 1][j];  



		alpha = A - B * g[j + 1];
		g[j] = C / alpha;
		v[j] = (F - B * v[j + 1]) / alpha;
	}

	//  Bottom boundary
	B = (j_minus_inviscid_Jacobians[i][1] + j_viscous_Jacobians[i][1]) * (grid.jArea(i, 1));

	A = grid.Volume(i, 0) / dt * identity()
		- (i_minus_inviscid_Jacobians[i][0] + i_viscous_Jacobians[i][0]) * grid.iArea(i, 0)
		+ (i_plus_inviscid_Jacobians[i + 1][0] - i_viscous_Jacobians[i + 1][0]) * grid.iArea(i + 1, 0)
		- (Ebi * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0] + j_viscous_Jacobians[i][0]) * grid.jArea(i, 0) 
		+ (j_plus_inviscid_Jacobians[i][1] - j_viscous_Jacobians[i][1]) * grid.jArea(i, 1); 

	F = i_Fluxes[i][0] * (grid.iArea(i, 0))
		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1))
		+ (i_plus_inviscid_Jacobians[i][0] - i_viscous_Jacobians[i][0]) * (grid.iArea(i, 0)) * dU_old[i - 1][0]
		+ (i_minus_inviscid_Jacobians[i + 1][0] + i_viscous_Jacobians[i + 1][0]) * (-grid.iArea(i + 1, 0)) * dU_old[i + 1][0]; 


	alpha = A - B * g[1];
	v[0] = (F - B * v[1]) / alpha;

	// Calculate dU
	dU_new[i][0] = v[0];

	for (int j = 1; j < Ny; ++j) {
		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];
	} 
}

void Solver::solve_right_line_viscous() {

	int i = Nx - 1;
	Matrix Eti, Ebi, Eri; 

	static Matrix A; // Relates to U(j+1) 
	static Matrix B; // Relates to U(j) 
	static Matrix C; // Relates to U(j-1) 
	static Vector F; // Right hand side  

	// Intermediate matrices
	static Matrix alpha;
	static array<array<double, 4>, Ny> v;
	static array<array<array<double, 4>, 4>, Ny> g;

	// Grab boundary values
	Ebi = inviscid_boundary_2D_E(BoundaryType.bottom, U[i][0], grid.jNorms(i, 0));
	Eti = inviscid_boundary_2D_E(BoundaryType.top, U[i][Ny - 1], grid.jNorms(i, Ny)); 
	Eri = inviscid_boundary_2D_E(BoundaryType.right, U[i][Ny - 1], grid.iNorms(i, Ny - 1)); 

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity()
		- (i_minus_inviscid_Jacobians[i][Ny - 1] + i_viscous_Jacobians[i][Ny - 1]) * grid.iArea(i, Ny - 1)																// Left
		+ (i_plus_inviscid_Jacobians[i + 1][Ny - 1] + Eri * i_minus_inviscid_Jacobians[i + 1][Ny - 1] + i_viscous_Jacobians[i + 1][Ny - 1]) * grid.iArea(i + 1, Ny - 1)	// Right
		- (j_minus_inviscid_Jacobians[i][Ny - 1] + j_viscous_Jacobians[i][Ny - 1]) * grid.jArea(i, Ny - 1)																// Bottom
		+ (j_plus_inviscid_Jacobians[i][Ny] + Eti * j_minus_inviscid_Jacobians[i][Ny] + j_viscous_Jacobians[i][Ny]) * grid.jArea(i, Ny);		 						// Top 

	C = (j_plus_inviscid_Jacobians[i][Ny - 1] - j_viscous_Jacobians[i][Ny - 1]) * (-grid.jArea(i, Ny - 1)); 

	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1))
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny))
		+ (i_plus_inviscid_Jacobians[i][Ny - 1] - i_viscous_Jacobians[i][Ny - 1]) * (grid.iArea(i, Ny - 1)) * dU_old[i - 1][Ny - 1];

	alpha = A;
	v[Ny - 1] = F / alpha;
	g[Ny - 1] = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		Eri = inviscid_boundary_2D_E(BoundaryType.right, U[i][j], grid.iNorms(i, j));	

		B = (j_minus_inviscid_Jacobians[i][j + 1] + j_viscous_Jacobians[i][j + 1]) * (grid.jArea(i, j + 1));

		A = grid.Volume(i, j) / dt * identity()
			- (i_minus_inviscid_Jacobians[i][j] + i_viscous_Jacobians[i][j]) * grid.iArea(i, j)															// Left
			+ (i_plus_inviscid_Jacobians[i + 1][j] + Eri * i_minus_inviscid_Jacobians[i + 1][j] + i_viscous_Jacobians[i + 1][j]) * grid.iArea(i + 1, j)	// Right
			- (j_minus_inviscid_Jacobians[i][j] + j_viscous_Jacobians[i][j]) * grid.jArea(i, j)															// Bottom
			+ (j_plus_inviscid_Jacobians[i][j + 1] - j_viscous_Jacobians[i][j + 1]) * grid.jArea(i, j + 1);												// Top

		C = (j_plus_inviscid_Jacobians[i][j] - j_viscous_Jacobians[i][j]) * (-grid.jArea(i, j)); 

		F = i_Fluxes[i][j] * (grid.iArea(i, j))
			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j))
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1))
			+ (i_plus_inviscid_Jacobians[i][j] - i_viscous_Jacobians[i][j]) * (grid.iArea(i, j)) * dU_old[i - 1][j]; 



		alpha = A - B * g[j + 1];
		g[j] = C / alpha;
		v[j] = (F - B * v[j + 1]) / alpha;
	}
	//  Bottom boundary
	Eri = inviscid_boundary_2D_E(BoundaryType.right, U[i][0], grid.iNorms(i, 0));

	B = (j_minus_inviscid_Jacobians[i][1] + j_viscous_Jacobians[i][1]) * (grid.jArea(i, 1)); 

	A = grid.Volume(i, 0) / dt * identity()
		- (i_minus_inviscid_Jacobians[i][0] + i_viscous_Jacobians[i][0]) * grid.iArea(i, 0)															// Left
		+ (i_plus_inviscid_Jacobians[i + 1][0] + Eri * i_minus_inviscid_Jacobians[i + 1][0] + i_viscous_Jacobians[i + 1][0]) * grid.iArea(i + 1, 0)	// Right 
		- (Ebi * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0] + j_viscous_Jacobians[i][0]) * grid.jArea(i, 0)					// Bottom 
		+ (j_plus_inviscid_Jacobians[i][1] - j_viscous_Jacobians[i][1]) * grid.jArea(i, 1);															// Top 

	F = i_Fluxes[i][0] * (grid.iArea(i, 0)) 
		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0)) 
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1))
		+ (i_plus_inviscid_Jacobians[i][0] - i_viscous_Jacobians[i][0]) * (grid.iArea(i, 0)) * dU_old[i - 1][0]; 


	alpha = A - B * g[1];
	v[0] = (F - B * v[1]) / alpha;

	// Calculate dU
	dU_new[i][0] = v[0];

	for (int j = 1; j < Ny; ++j) {
		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];
	}

}


void Solver::compute_inner_residual() { 
	
	inner_residual = 0.0;
	double res = 0.0; 

	static Vector A; // Relates to U(j+1) 
	static Vector B; // Relates to U(j) 
	static Vector C; // Relates to U(j-1) 
	static double F; // Right hand side  

	Vector one = { 1, 0, 0, 0 };

	// Middle Boundary
	for (int i = 1; i < Nx - 1; ++i) {
		for (int j = 1; j < Ny - 1; ++j) {

			B = j_minus_inviscid_Jacobians[i][j + 1][0] * (grid.jArea(i, j + 1)); 

			A = grid.Volume(i, j) / dt * one 
				- i_minus_inviscid_Jacobians[i][j][0] * grid.iArea(i, j) 
				+ i_plus_inviscid_Jacobians[i + 1][j][0] * grid.iArea(i + 1, j)
				- j_minus_inviscid_Jacobians[i][j][0] * grid.jArea(i, j) 
				+ j_plus_inviscid_Jacobians[i][j + 1][0] * grid.jArea(i, j + 1);

			C = j_plus_inviscid_Jacobians[i][j][0] * (-grid.jArea(i, j)); 

			F = i_Fluxes[i][j][0] * (grid.iArea(i, j)) 
				+ i_Fluxes[i + 1][j][0] * (-grid.iArea(i + 1, j)) 
				+ j_Fluxes[i][j][0] * (grid.jArea(i, j)) 
				+ j_Fluxes[i][j + 1][0] * (-grid.jArea(i, j + 1)) 
				+ i_plus_inviscid_Jacobians[i][j][0] * (grid.iArea(i, j)) * dU_old[i - 1][j]
				+ i_minus_inviscid_Jacobians[i + 1][j][0] * (-grid.iArea(i + 1, j)) * dU_old[i + 1][j];   

			res = B * dU_new[i][j + 1] + A * dU_new[i][j] + C * dU_new[i][j - 1] - F;
			inner_residual += res * res;  

		}
	}
	inner_residual = sqrt(inner_residual); 
}

void Solver::compute_outer_residual() {  

	outer_residual = 0.0;  
	double intres;

	for (int i = 1; i < Nx - 1; ++i) {
		for (int j = 1; j < Ny - 1; ++j) {
			intres = (i_Fluxes[i][j][0] * (-grid.iArea(i, j))			// Left Side		   
				+ i_Fluxes[i + 1][j][0] * (grid.iArea(i + 1, j))		// Right Side 
				+ j_Fluxes[i][j][0] * (-grid.jArea(i, j))			// Bottom Side  
				+ j_Fluxes[i][j + 1][0] * (grid.jArea(i, j + 1)))/ grid.Volume(i, j); 

			outer_residual = outer_residual + intres * intres; 
		
		}
	}

	outer_residual = sqrt(outer_residual);
}

void Solver::write_2d_csv(const string& filename) {

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

void Solver::write_1d_csv(const string& filename) {
	double a;

	array<double, Ny> density, u_velocity, v_velocity, pressure, Mach, Temperature; 
	Vector Primitives;


	for (int j = 0; j < Ny; ++j) {

		Primitives = constoPrim(U[Nx/2][j]);  
		density[j] = Primitives[0]; 
		pressure[j] = Primitives[3]; 
		a = sqrt(gamma * pressure[j] / density[j]); 
		u_velocity[j] = Primitives[1]; 
		v_velocity[j] = Primitives[2]; 
		Mach[j] = sqrt(u_velocity[j] * u_velocity[j] + v_velocity[j] * v_velocity[j]) / a; 
		Temperature[j] = pressure[j] / (density[j] * R); 
	}
	
	ofstream file(filename);

	file << "density, u velocity, v velocity, pressure, Mach, Temperature, y_points" << endl;

	for (int j = 0; j < Ny; ++j) { 
		file << density[j]/INLET.rho << ", " << u_velocity[j]/INLET.u << ", " << v_velocity[j]/INLET.v << ", " << pressure[j]/INLET.p 
			<< ", " << Mach[j]/INLET.M << ", " << Temperature[j]/INLET.T << ", "  << grid.Center(Nx / 2, j).y << endl;  
	}
	

	file.close();
	cout << "1D File saved successfully as \"" << filename << "\"" << endl; 

}

