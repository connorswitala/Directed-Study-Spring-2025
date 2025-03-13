// 2D FVS Library.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "framework.h"
#include "2DFVSLibrary.h"


Vector primtoCons(const Vector& V) {

	Vector U(4);
	U[0] = V[0];
	U[1] = V[0] * V[1];
	U[2] = V[0] * V[2];
	U[3] = V[3] / (gamma - 1) + 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2]);

	return U;
}

Vector constoPrim(const Vector& U) {

	static Vector V(4);
	V[0] = U[0];
	V[1] = U[1] / U[0];
	V[2] = U[2] / U[0];
	V[3] = (U[3] - 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2])) * (gamma - 1);

	return V;
}

Solver::Solver(const int Nx, const int Ny, const inlet_conditions& INLET, Grid& grid, BoundaryConditions BoundaryType, double CFL, double Tw, int& progress_update) : Nx(Nx), Ny(Ny), 
Tw(Tw), INLET(INLET), V_inlet(V_inlet), U_inlet(U_inlet), grid(grid), BoundaryType(BoundaryType), U(U), dU_new(dU_new), dU_old(dU_old), i_Fluxes(i_Fluxes), j_Fluxes(j_Fluxes), 
i_plus_inviscid_Jacobians(i_plus_inviscid_Jacobians), j_plus_inviscid_Jacobians(j_plus_inviscid_Jacobians), i_minus_inviscid_Jacobians(i_minus_inviscid_Jacobians), 
j_minus_inviscid_Jacobians(j_minus_inviscid_Jacobians), i_viscous_Jacobians(i_viscous_Jacobians), j_viscous_Jacobians(j_viscous_Jacobians),  gridtype(gridtype), dt(dt), 
CFL(CFL), inner_residual(inner_residual), outer_residual(outer_residual), progress_update(progress_update){

	V_inlet = Vector(4); U_inlet = Vector(4); 
	U = Tensor(Nx + 4, Matrix(Ny + 4, Vector(4, 0.0)));  
	dU_new = Tensor(Nx, Matrix(Ny, Vector(4, 0.0)));  
	dU_old = Tensor(Nx, Matrix(Ny, Vector(4, 0.0)));   

	i_Fluxes = Tensor(Nx + 1, Matrix(Ny, Vector(4, 0.0)));
	j_Fluxes = Tensor(Nx, Matrix(Ny + 1, Vector(4, 0.0)));

	i_plus_inviscid_Jacobians = Tesseract(Nx + 1, Tensor(Ny, Matrix(4, Vector(4, 0.0))));
	i_minus_inviscid_Jacobians = Tesseract(Nx + 1, Tensor(Ny, Matrix(4, Vector(4, 0.0)))); 
	i_viscous_Jacobians = Tesseract(Nx + 1, Tensor(Ny, Matrix(4, Vector(4, 0.0)))); 

	j_plus_inviscid_Jacobians = Tesseract(Nx, Tensor(Ny + 1, Matrix(4, Vector(4, 0.0))));
	j_minus_inviscid_Jacobians = Tesseract(Nx, Tensor(Ny + 1, Matrix(4, Vector(4, 0.0))));
	j_viscous_Jacobians = Tesseract(Nx, Tensor(Ny + 1, Matrix(4, Vector(4, 0.0))));  
	

	outer_residual = 1.0;
	inner_residual = 1.0; 

	V_inlet = { INLET.rho, INLET.u, INLET.v, INLET.p };

	gridtype;   
	if (dynamic_cast<RampGrid*>(&grid)) gridtype = "Ramp";  
	else if (dynamic_cast<CylinderGrid*>(&grid)) gridtype = "Cylinder";  
	else if (dynamic_cast<FlatPlateGrid*>(&grid)) gridtype = "Flat Plate"; 
	else if (dynamic_cast<DoubleConeGrid*>(&grid)) gridtype = "Double Cone"; 
	else if (dynamic_cast<MirroredGrid*>(&grid)) gridtype = "Mirrored Double Ramp";
	else gridtype = "Unknown";  


	U_inlet = primtoCons(V_inlet); 

	for (int i = 0; i < Nx + 4; ++i) {
		for (int j = 0; j < Ny + 4; ++j) {
			U[i][j] = U_inlet;  
		}
	}

};   

Vector Solver::constoViscPrim(const Vector& U) { 
	Vector result(4);
	result[0] = U[0];
	result[1] = U[1] / U[0];
	result[2] = U[2] / U[0];
	result[3] = computeTemperature(U); 
	return result;  
}

Matrix Solver::inviscid_boundary_2D_E(BoundaryCondition type, const Vector& U, const Point& normals) {     

	Matrix E = zeros(4, 4);  	

	switch (type) {

	case BoundaryCondition::Inlet: 
		return E;  

	case BoundaryCondition::Outlet: 
		
		return identity(4);  

	case BoundaryCondition::IsothermalWall:  
	 
		E = { {1, 0, 0, 0}, 
			  {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
			  { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
			  { 0, 0, 0, 1 }
		};

		return E;

	case BoundaryCondition::AdiabaticWall:

		E = { {1, 0, 0, 0}, 
			  {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
			  { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
			  { 0, 0, 0, 1 }
		};

		return E;

	case BoundaryCondition::Symmetry: 		

		E = { {1, 0, 0, 0},
			  {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
			  { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
			  { 0, 0, 0, 1 }
		};

		return E; 
	default:
		throw invalid_argument("Unknown boundary condition type.");
	}

}

Vector Solver::inviscid_boundary_2D_U(BoundaryCondition type, const Vector& U, const Point& normals) {

	Vector ghost(4);
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

	Matrix E = zeros(4, 4); 
	double Ti; 

	switch (type) {

	case BoundaryCondition::Inlet:
		return E; 

	case BoundaryCondition::Outlet:

		return identity(4);

	case BoundaryCondition::IsothermalWall:

		Ti = computeTemperature(U); 

		E = { {2 * Ti/Tw - 1 , 0, 0, 2 * U[0]/Tw}, 
			  { 0, -1, 0, 0 },
			  { 0, 0, -1, 0},
			  { 0, 0, 0, -1 }
		};

		return E;

	case BoundaryCondition::AdiabaticWall:

		Ti = computeTemperature(U);

		E = { {2 * Ti / Tw - 1 , 0, 0, 2 * U[0] / Tw}, 
			  { 0, -1, 0, 0 },
			  { 0, 0, -1, 0},
			  { 0, 0, 0, 1 }
		};

		return E; 

	case BoundaryCondition::Symmetry:

		E = { {1, 0, 0, 0},
			  {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
			  { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
			  { 0, 0, 0, 1 }
		};

		return E;
	default:
		throw invalid_argument("Unknown boundary condition type.");
	}

}

Vector Solver::viscous_boundary_2D_U(BoundaryCondition type, const Vector& U, const Point& normals) {

	Vector ghost(4);
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

	Vector V(4); 
	double dx, dy, c, dt_old; 
	dt = 1; 

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			V = constoPrim(U[i + 2][j + 2]); 
			dx = min(grid.jArea(i, j), grid.jArea(i, j + 1));
			dy = min(grid.iArea(i, j), grid.iArea(i + 1, j));
			c = sqrt(gamma * V[3] / V[0]);

			dt_old = CFL / (fabs(V[1] / dx) + fabs(V[2] / dy) + c * sqrt(1 / (dx * dx) + 1 / (dy * dy))); 

			if (dt_old < dt) dt = dt_old; 

		}
	}
}

void Solver::solve_inviscid () { 

	cout << "\033[33mRunning Inviscid DPLR for " << Nx << " by " << Ny << " " << gridtype << " with a CFL of " << fixed << setprecision(2) << CFL << "...\033[0m" << "\n\n";
	string filename = "Inviscid " + to_string(Nx) + "x" + to_string(Ny) + "_" + gridtype + "_Solution.csv"; 
	
	auto start = TIME; 
	int counter = 0;

	Chemistry chem(INLET); 

	while (outer_residual >= 1e-6) { 

		compute_dt();

		#pragma omp parallel for
		for (int i = 0; i < Nx; ++i) {
			for (int j = 0; j < Ny; ++j) {
				U[i + 2][j + 2] = chem.compute_equilibrium(U[i + 2][j + 2]);  
			}
		}
	 
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

		if (counter % 100 == 0) { 
			write_2d_csv(filename); 
		}


	}

	write_2d_csv(filename);  

	auto end = TIME;
	DURATION duration = end - start; 
	cout << "Program completed in " << duration.count() << " seconds." << endl;


}  

void Solver::solve_viscous() {

	cout << "\033[33mRunning Viscous DPLR for " << Nx << " by " << Ny << " " << gridtype << " with a CFL of " << fixed << setprecision(2) << CFL << "...\033[0m" << "\n\n"; 
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
	cout << "Program complete in " << duration.count() << " seconds." << endl;



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

		fill_ghost_cells_inviscid(); 

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
				U[i + 2][j + 2][k] += dU_old[i][j][k];     
			}
		}
	}

}

void Solver::fill_ghost_cells_inviscid() { 

	for (int j = 0; j < Ny; ++j) { 
		U[0][j + 2] = inviscid_boundary_2D_U(BoundaryType.left, U[2][j + 2], grid.iNorms(0, j));
		U[1][j + 2] = inviscid_boundary_2D_U(BoundaryType.left, U[2][j + 2], grid.iNorms(0, j));
		U[Nx + 2][j + 2] = inviscid_boundary_2D_U(BoundaryType.right, U[Nx + 1][j + 2], grid.iNorms(Nx, j));    
		U[Nx + 3][j + 2] = inviscid_boundary_2D_U(BoundaryType.right, U[Nx + 1][j + 2], grid.iNorms(Nx, j));   
	}

	for (int i = 0; i < Nx; ++i) {
		U[i + 2][0] = inviscid_boundary_2D_U(BoundaryType.bottom, U[i + 2][2], grid.jNorms(i, 0));    
		U[i + 2][1] = inviscid_boundary_2D_U(BoundaryType.bottom, U[i + 2][2], grid.jNorms(i, 0));     
		U[i + 2][Ny + 2] = inviscid_boundary_2D_U(BoundaryType.top, U[i + 2][Ny + 1], grid.jNorms(i, Ny));     
		U[i + 2][Ny + 3] = inviscid_boundary_2D_U(BoundaryType.top, U[i + 2][Ny + 1], grid.jNorms(i, Ny));      
	}
}

void Solver::compute_inviscid_jacobians() {

	static Vector Ui(4), Uii(4), Up(4), Um(4), V1_Plus(4), V2_Plus(4), V1_Minus(4), V2_Minus(4), n(4), m(4); 
	Inviscid_State Si, Sii; 
	double g = 5.72;
	double weight, dp, pi, pii, pe, nx, ny, l, lc, lt;
	pe = (gamma - 1);  
	Matrix X(4, Vector(4)), Y(4, Vector(4)); 

	// Calculate Jacobians and Explicit fluxes for i-faces on left boundary 
	for (int j = 0; j < Ny; ++j) {  
				
		Uii = U[2][j + 2];    
		Ui = U[1][j + 2];    

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

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity(4);

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


		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity(4);   
		i_plus_inviscid_Jacobians[0][j] = X;  
		i_minus_inviscid_Jacobians[0][j] = Y;  

		//Ui = Ui + 0.5 * minmod(U[1][j + 2] - U[0][j + 2], U[2][j + 2] - U[1][j + 2]); 
		//Uii = Uii - 0.5 * minmod(U[2][j + 2] - U[1][j + 2], U[3][j + 2] - U[2][j + 2]); 

		i_Fluxes[0][j] = X * Ui + Y * Uii;     
	}

	// Calculate Jacobians and Explicit fluxes for i-faces on right boundary
	for (int j = 0; j < Ny; ++j) {

		Ui = U[Nx + 1][j + 2]; 
		Uii = U[Nx + 2][j + 2]; 

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

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity(4);

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

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity(4);
		

		i_plus_inviscid_Jacobians[Nx][j] = X; 
		i_minus_inviscid_Jacobians[Nx][j] = Y; 

		//Ui = Ui + 0.5 * minmod(U[Nx + 1][j + 2] - U[Nx][j + 2], U[Nx + 2][j + 2] - U[Nx + 1][j + 2]);
		//Uii = Uii - 0.5 * minmod(U[Nx + 2][j + 2] - U[Nx + 1][j + 2], U[Nx + 3][j + 2] - U[Nx + 2][j + 2]);

		i_Fluxes[Nx][j] = X * Ui + Y * Uii;   
	}

	// Calculate Jacobians and Explicit fluxes for j-faces on bottom boundary
	for (int i = 0; i < Nx; ++i) { 

		Uii = U[i + 2][2]; 
		Ui = U[i + 2][1];  

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

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity(4);

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

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity(4);

		j_plus_inviscid_Jacobians[i][0] = X;
		j_minus_inviscid_Jacobians[i][0] = Y;   

		//Ui = Ui + 0.5 * minmod(U[i + 2][1] - U[i + 2][0], U[i + 2][2] - U[i + 2][1]);
		//Uii = Uii - 0.5 * minmod(U[i + 2][2] - U[i + 2][1], U[i + 2][3] - U[i + 2][2]);

		j_Fluxes[i][0] = X * Ui + Y * Uii;   

	}

	// Calculate Jacobians and Explicit fluxes for j-faces on top boundary
	for (int i = 0; i < Nx; ++i) {

		Ui = U[i + 2][Ny + 1]; 
		Uii = U[i + 2][Ny + 2]; 

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

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity(4);

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

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity(4); 

		j_plus_inviscid_Jacobians[i][Ny] = X; 
		j_minus_inviscid_Jacobians[i][Ny] = Y; 

		//Ui = Ui + 0.5 * minmod(U[i + 2][Ny + 1] - U[i + 2][Ny], U[i + 2][Ny + 2] - U[i + 2][Ny + 1]);
		//Uii = Uii - 0.5 * minmod(U[i + 2][Ny + 2] - U[i + 2][Ny + 1], U[i + 2][Ny + 3] - U[i + 2][Ny + 2]);

		j_Fluxes[i][Ny] = X * Ui + Y * Uii;   

	}

	// Inner i-faces  
	for (int i = 1; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {	

			Ui = U[i + 1][j + 2]; 
			Uii = U[i + 2][j + 2]; 

			nx = grid.iNorms(i, j).x;
			ny = grid.iNorms(i, j).y;

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

			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity(4);
			

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
			
			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity(4); 
			 
			i_plus_inviscid_Jacobians[i][j] = X;
			i_minus_inviscid_Jacobians[i][j] = Y;

			//Ui = Ui + 0.5 * minmod(U[i + 1][j + 2] - U[i][j + 2], U[i + 2][j + 2] - U[i + 1][j + 2]);    
			//Uii = Uii - 0.5 * minmod(U[i + 2][j + 2] - U[i + 1][j + 2], U[i + 3][j + 2] - U[i + 2][j + 2]); 

			i_Fluxes[i][j] = X * Ui + Y * Uii;   
		}
	}

	// Inner j-faces 
	for (int i = 0; i < Nx; ++i) {
		for (int j = 1; j < Ny; ++j) {
			
			nx = grid.jNorms(i, j).x;
			ny = grid.jNorms(i, j).y;

			Ui = U[i + 2][j + 1];  
			Uii = U[i + 2][j + 2];   
		
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

			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity(4);

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

			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity(4); 

			j_plus_inviscid_Jacobians[i][j] = X; 
			j_minus_inviscid_Jacobians[i][j] = Y;

			//Ui = Ui + 0.5 * minmod(U[i + 2][j + 1] - U[i + 2][j], U[i + 2][j + 2] - U[i + 2][j + 1]);  
			//Uii = Uii - 0.5 * minmod(U[i + 2][j + 2] - U[i + 2][j + 1], U[i + 2][j + 3] - U[i + 2][j + 2]); 

			j_Fluxes[i][j] = X * Ui + Y * Uii;   
		}
	}

} 

void Solver::compute_viscous_jacobians() {  

	static Vector Ui(4), Uii(4), Up(4), Um(4), V1_Plus(4), V2_Plus(4), V1_Minus(4), V2_Minus(4), n(4), m(4), dv(4), W(4);
	Viscous_State Si, Sii; 
	double g = 5.72;
	double weight, dp, pi, pii, pe, nx, ny, lp, lcp, ltp, lm, lcm, ltm, rho, u, v, T, mu, lambda, k, dn, dx, dy;  
	pe = (gamma - 1);
	Matrix M(4, Vector(4)), N(4, Vector(4)), X(4, Vector(4)), Y(4, Vector(4)), Elv(4, Vector(4)), Erv(4, Vector(4)), Ebv(4, Vector(4)), Etv(4, Vector(4));

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

		M[0] = { 0, 0, 0, 0 };
		M[1][0] = 0;
		M[1][1] = (lambda + 2 * mu) * nx * nx + mu * ny * ny; 
		M[1][2] = (lambda + mu) * nx * ny;
		M[1][3] = 0;
		M[2][0] = 0;
		M[2][1] = (mu + lambda) * nx * ny;
		M[2][2] = mu * nx * nx + (lambda + 2 * mu) * ny * ny;
		M[2][3] = 0;
		M[3][0] = 0;
		M[3][1] = u * (lambda + 2 * mu) * nx * nx + (v * mu + v * lambda) * nx * ny + u * mu * ny * ny;  
		M[3][2] = v * mu * nx * nx + (u * lambda + u * mu) * nx * ny + v * (lambda + 2 * mu) * ny * ny; 
		M[3][3] = k * (nx * nx + ny * ny);  


		N = {
			{ {1, 0, 0, 0},
			{-u / rho, 1 / rho, 0, 0},
			{-v / rho, 0, 1 / rho, 0},
			{ (u * u + v * v) / (2 * cv * rho) - T / rho, -u / (cv * rho), -v / (cv * rho), 1 / (cv * rho)}
			}
		};
	
		dx = grid.Center(1, j).x - grid.Center(0, j).x;
		dy = grid.Center(1, j).y - grid.Center(0, j).y; 
		dn = sqrt(dx * dx + dy * dy); 

		dv = constoViscPrim(Uii) - constoViscPrim(Ui);

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

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(4);

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


		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(4);

		i_viscous_Jacobians[0][j] = M/(-dn) * (identity(4) - Elv) * N;  

		//displayMatrix(Z);  
		i_plus_inviscid_Jacobians[0][j] = X;
		i_minus_inviscid_Jacobians[0][j] = Y;
		i_Fluxes[0][j] = X * Ui + Y * Uii - M * dv / dn;  
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

		M[0] = { 0, 0, 0, 0 };
		M[1][0] = 0;
		M[1][1] = (lambda + 2 * mu) * nx * nx + mu * ny * ny;
		M[1][2] = (lambda + mu) * nx * ny;
		M[1][3] = 0;
		M[2][0] = 0;
		M[2][1] = (mu + lambda) * nx * ny;
		M[2][2] = mu * nx * nx + (lambda + 2 * mu) * ny * ny;
		M[2][3] = 0;
		M[3][0] = 0;
		M[3][1] = u * (lambda + 2 * mu) * nx * nx + (v * mu + v * lambda) * nx * ny + u * mu * ny * ny;
		M[3][2] = v * mu * nx * nx + (u * lambda + u * mu) * nx * ny + v * (lambda + 2 * mu) * ny * ny;
		M[3][3] = k * (nx * nx + ny * ny);

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

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(4);

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

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(4);

		i_viscous_Jacobians[Nx][j] = M/(-dn) * (Erv - identity(4)) * N;  
		i_plus_inviscid_Jacobians[Nx][j] = X;
		i_minus_inviscid_Jacobians[Nx][j] = Y;
		i_Fluxes[Nx][j] = X * Ui + Y * Uii - M * dv / dn; 
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

		M[0] = { 0, 0, 0, 0 };
		M[1][0] = 0;
		M[1][1] = (lambda + 2 * mu) * nx * nx + mu * ny * ny;
		M[1][2] = (lambda + mu) * nx * ny;
		M[1][3] = 0;
		M[2][0] = 0;
		M[2][1] = (mu + lambda) * nx * ny;
		M[2][2] = mu * nx * nx + (lambda + 2 * mu) * ny * ny;
		M[2][3] = 0;
		M[3][0] = 0;
		M[3][1] = u * (lambda + 2 * mu) * nx * nx + (v * mu + v * lambda) * nx * ny + u * mu * ny * ny;
		M[3][2] = v * mu * nx * nx + (u * lambda + u * mu) * nx * ny + v * (lambda + 2 * mu) * ny * ny;
		M[3][3] = k * (nx * nx + ny * ny);

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

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(4);

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

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(4);

		j_viscous_Jacobians[i][0] = M/(-dn) * (identity(4) - Ebv) * N;   
		j_plus_inviscid_Jacobians[i][0] = X;
		j_minus_inviscid_Jacobians[i][0] = Y;
		j_Fluxes[i][0] = X * Ui + Y * Uii - M * dv / dn;

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

		M[0] = { 0, 0, 0, 0 };
		M[1][0] = 0;
		M[1][1] = (lambda + 2 * mu) * nx * nx + mu * ny * ny;
		M[1][2] = (lambda + mu) * nx * ny;
		M[1][3] = 0;
		M[2][0] = 0;
		M[2][1] = (mu + lambda) * nx * ny;
		M[2][2] = mu * nx * nx + (lambda + 2 * mu) * ny * ny;
		M[2][3] = 0;
		M[3][0] = 0;
		M[3][1] = u * (lambda + 2 * mu) * nx * nx + (v * mu + v * lambda) * nx * ny + u * mu * ny * ny;
		M[3][2] = v * mu * nx * nx + (u * lambda + u * mu) * nx * ny + v * (lambda + 2 * mu) * ny * ny;
		M[3][3] = k * (nx * nx + ny * ny);

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

		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(4);

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

		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(4);

		j_viscous_Jacobians[i][Ny] = M/(-dn) * (Etv - identity(4)) * N; 
		j_plus_inviscid_Jacobians[i][Ny] = X;
		j_minus_inviscid_Jacobians[i][Ny] = Y;
		j_Fluxes[i][Ny] = X * Ui + Y * Uii - M * dv / dn;

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

			M[0] = { 0, 0, 0, 0 };
			M[1][0] = 0;
			M[1][1] = (lambda + 2 * mu) * nx * nx + mu * ny * ny;
			M[1][2] = (lambda + mu) * nx * ny;
			M[1][3] = 0;
			M[2][0] = 0;
			M[2][1] = (mu + lambda) * nx * ny;
			M[2][2] = mu * nx * nx + (lambda + 2 * mu) * ny * ny;
			M[2][3] = 0;
			M[3][0] = 0;
			M[3][1] = u * (lambda + 2 * mu) * nx * nx + (v * mu + v * lambda) * nx * ny + u * mu * ny * ny;
			M[3][2] = v * mu * nx * nx + (u * lambda + u * mu) * nx * ny + v * (lambda + 2 * mu) * ny * ny;
			M[3][3] = k * (nx * nx + ny * ny);

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

			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(4);


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

			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(4);

			i_viscous_Jacobians[i][j] = M/(-dn) * N;  
			i_plus_inviscid_Jacobians[i][j] = X;
			i_minus_inviscid_Jacobians[i][j] = Y;
			i_Fluxes[i][j] = X * Ui + Y * Uii - M * dv / dn; 
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

			M[0] = { 0, 0, 0, 0 };
			M[1][0] = 0;
			M[1][1] = (lambda + 2 * mu) * nx * nx + mu * ny * ny;
			M[1][2] = (lambda + mu) * nx * ny;
			M[1][3] = 0;
			M[2][0] = 0;
			M[2][1] = (mu + lambda) * nx * ny;
			M[2][2] = mu * nx * nx + (lambda + 2 * mu) * ny * ny;
			M[2][3] = 0;
			M[3][0] = 0;
			M[3][1] = u * (lambda + 2 * mu) * nx * nx + (v * mu + v * lambda) * nx * ny + u * mu * ny * ny;
			M[3][2] = v * mu * nx * nx + (u * lambda + u * mu) * nx * ny + v * (lambda + 2 * mu) * ny * ny;
			M[3][3] = k * (nx * nx + ny * ny);


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

			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(4);

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

			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(4);

			j_viscous_Jacobians[i][j] = M/(-dn) * N; 
			j_plus_inviscid_Jacobians[i][j] = X;
			j_minus_inviscid_Jacobians[i][j] = Y;
			j_Fluxes[i][j] = X * Ui + Y * Uii - M * dv/dn;  
		}
	}

}

void Solver::solve_left_line_inviscid() {

	int i = 0; 
	Matrix Et(4, Vector(4)), Eb(4, Vector(4)), El(4, Vector(4));
	static Matrix A(4, Vector(4)); // Relates to U(j+1) 
	static Matrix B(4, Vector(4)); // Relates to U(j) 
	static Matrix C(4, Vector(4)); // Relates to U(j-1) 
	static Vector F(4); // Right hand side  

	// Intermediate matrices
	static Matrix alpha(4, Vector(4));
	static Matrix v(Ny, Vector(4));
	static Tensor g(Ny, Matrix(4, Vector(4))); 

	// Grab boundary values
	Eb = inviscid_boundary_2D_E(BoundaryType.bottom, U[i + 2][2], grid.jNorms(i, 0));      
	Et = inviscid_boundary_2D_E(BoundaryType.top, U[i + 2][Ny + 1], grid.jNorms(i, Ny));     
	El = inviscid_boundary_2D_E(BoundaryType.left, U[i + 2][Ny + 1], grid.iNorms(i, Ny - 1)); 

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity(4)
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

		El = inviscid_boundary_2D_E(BoundaryType.left, U[i + 2][j + 2], grid.iNorms(i, j));   

		B = j_minus_inviscid_Jacobians[i][j + 1] * (grid.jArea(i, j + 1));

		A = grid.Volume(i, j) / dt * identity(4)
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

	El = inviscid_boundary_2D_E(BoundaryType.left, U[i + 2][2], grid.iNorms(i, 0));  

	B = j_minus_inviscid_Jacobians[i][1] * (grid.jArea(i, 1));

	A = grid.Volume(i, 0) / dt * identity(4)
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
	Matrix Et(4, Vector(4)), Eb(4, Vector(4));
	static Matrix A(4, Vector(4)); // Relates to U(j+1) 
	static Matrix B(4, Vector(4)); // Relates to U(j) 
	static Matrix C(4, Vector(4)); // Relates to U(j-1) 
	static Vector F(4); // Right hand side  

	// Intermediate matrices
	static Matrix alpha(4, Vector(4));
	static Matrix v(Ny, Vector(4));
	static Tensor g(Ny, Matrix(4, Vector(4)));

	// Grab boundary values
	Eb = inviscid_boundary_2D_E(BoundaryType.bottom, U[i + 2][2], grid.jNorms(i, 0));  
	Et = inviscid_boundary_2D_E(BoundaryType.top, U[i + 2][Ny + 1], grid.jNorms(i, Ny));    

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity(4) 
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

		A = grid.Volume(i, j) / dt * identity(4) 
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

	A = grid.Volume(i, 0) / dt * identity(4)
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
	Matrix Et(4, Vector(4)), Eb(4, Vector(4)), Er(4, Vector(4));

	static Matrix A(4, Vector(4)); // Relates to U(j+1) 
	static Matrix B(4, Vector(4)); // Relates to U(j) 
	static Matrix C(4, Vector(4)); // Relates to U(j-1) 
	static Vector F(4); // Right hand side  

	// Intermediate matrices
	static Matrix alpha(4, Vector(4));
	static Matrix v(Ny, Vector(4));
	static Tensor g(Ny, Matrix(4, Vector(4)));

	// Grab boundary values
	Eb = inviscid_boundary_2D_E(BoundaryType.bottom, U[i + 2][2], grid.jNorms(i, 0)); 
	Et = inviscid_boundary_2D_E(BoundaryType.top, U[i + 2][Ny + 1], grid.jNorms(i, Ny)); 
	Er = inviscid_boundary_2D_E(BoundaryType.right, U[i + 2][Ny + 1], grid.iNorms(i + 1, Ny - 1));  

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity(4)
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

		Er = inviscid_boundary_2D_E(BoundaryType.right, U[i + 2][j + 2], grid.iNorms(i + 1, j));   

		B = j_minus_inviscid_Jacobians[i][j + 1] * (grid.jArea(i, j + 1));

		A = grid.Volume(i, j) / dt * identity(4)
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
	Er = inviscid_boundary_2D_E(BoundaryType.right, U[i + 2][2], grid.iNorms(i + 1, 0));  

	B = j_minus_inviscid_Jacobians[i][1] * (grid.jArea(i, 1));

	A = grid.Volume(i, 0) / dt * identity(4)
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
	Matrix Eti(4, Vector(4)), Ebi(4, Vector(4)), Eli(4, Vector(4));
	static Matrix A(4, Vector(4)); // Relates to U(j+1) 
	static Matrix B(4, Vector(4)); // Relates to U(j) 
	static Matrix C(4, Vector(4)); // Relates to U(j-1) 
	static Vector F(4); // Right hand side  

	// Intermediate matrices
	static Matrix alpha(4, Vector(4));
	static Matrix v(Ny, Vector(4));
	static Tensor g(Ny, Matrix(4, Vector(4)));

	// Grab boundary values
	Ebi = inviscid_boundary_2D_E(BoundaryType.bottom, U[i][0], grid.jNorms(i, 0));
	Eti = inviscid_boundary_2D_E(BoundaryType.top, U[i][Ny - 1], grid.jNorms(i, Ny));

	Eli = inviscid_boundary_2D_E(BoundaryType.left, U[i][Ny - 1], grid.iNorms(i, Ny - 1));

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity(4)
		- (i_minus_inviscid_Jacobians[i][Ny - 1] + Eli * i_plus_inviscid_Jacobians[i][Ny - 1] + i_viscous_Jacobians[i][Ny - 1]) * grid.iArea(i, Ny - 1)		// Left     
		+ (i_plus_inviscid_Jacobians[i + 1][Ny - 1] - i_viscous_Jacobians[i + 1][Ny - 1]) * grid.iArea(i + 1, Ny - 1)										// Right
		- (j_minus_inviscid_Jacobians[i][Ny - 1] + j_viscous_Jacobians[i][Ny - 1]) * grid.jArea(i, Ny - 1)		 											// Bottom
		+ (j_plus_inviscid_Jacobians[i][Ny] + Eti * j_minus_inviscid_Jacobians[i][Ny] + j_viscous_Jacobians[i][Ny]) * grid.jArea(i, Ny);					// Top  

	C = (j_plus_inviscid_Jacobians[i][Ny - 1] - j_viscous_Jacobians[i][Ny - 1]) * (-grid.jArea(i, Ny - 1)); 

	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
		+ (i_Fluxes[i + 1][Ny - 1] + (i_minus_inviscid_Jacobians[i + 1][Ny - 1] + i_viscous_Jacobians[i + 1][Ny - 1]) * dU_old[i + 1][Ny - 1]) * (-grid.iArea(i + 1, Ny - 1))
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny));

	alpha = A;
	v[Ny - 1] = F / alpha;
	g[Ny - 1] = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		Eli = inviscid_boundary_2D_E(BoundaryType.left, U[i][j], grid.iNorms(i, j)); 

		B = (j_minus_inviscid_Jacobians[i][j + 1] + j_viscous_Jacobians[i][j + 1]) * (grid.jArea(i, j + 1));

		A = grid.Volume(i, j) / dt * identity(4)
			- (i_minus_inviscid_Jacobians[i][j] + Eli * i_plus_inviscid_Jacobians[i][j] + i_viscous_Jacobians[i][j]) * grid.iArea(i, j)			// Left  
			+ (i_plus_inviscid_Jacobians[i + 1][j] - i_viscous_Jacobians[i + 1][j] ) * grid.iArea(i + 1, j)										// Right
			- (j_minus_inviscid_Jacobians[i][j] + j_viscous_Jacobians[i][j]) * grid.jArea(i, j)													// Bottom
			+ (j_plus_inviscid_Jacobians[i][j + 1] - j_viscous_Jacobians[i][j + 1]) * grid.jArea(i, j + 1);										// Top

		C = (j_plus_inviscid_Jacobians[i][j] - j_viscous_Jacobians[i][j]) * (-grid.jArea(i, j));  

		F = i_Fluxes[i][j] * (grid.iArea(i, j))
			+ (i_Fluxes[i + 1][j] + (i_minus_inviscid_Jacobians[i + 1][j] + i_viscous_Jacobians[i + 1][j]) * dU_old[i + 1][j]) * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j))
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1));



		alpha = A - B * g[j + 1];
		g[j] = C / alpha;
		v[j] = (F - B * v[j + 1]) / alpha;
	}

	//  Bottom boundary

	Eli = inviscid_boundary_2D_E(BoundaryType.left, U[i][0], grid.iNorms(i, 0));

	B = (j_minus_inviscid_Jacobians[i][1] + j_viscous_Jacobians[i][1]) * (grid.jArea(i, 1)); 

	A = grid.Volume(i, 0) / dt * identity(4)
		- (i_minus_inviscid_Jacobians[i][0] + Eli * i_plus_inviscid_Jacobians[i][0] + i_viscous_Jacobians[i][0]) * grid.iArea(i, 0)		// Left
		+ (i_plus_inviscid_Jacobians[i + 1][0] - i_viscous_Jacobians[i + 1][0]) * grid.iArea(i + 1, 0)									// Right 
		- (Ebi * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0] + j_viscous_Jacobians[i][0]) * grid.jArea(i, 0)		// Bottom
		+ (j_plus_inviscid_Jacobians[i][1] - j_viscous_Jacobians[i][1]) * grid.jArea(i, 1);												// Top

	F = i_Fluxes[i][0] * (grid.iArea(i, 0))
		+ (i_Fluxes[i + 1][0] + (i_minus_inviscid_Jacobians[i + 1][0] + i_viscous_Jacobians[i + 1][0]) * dU_old[i + 1][0]) * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1));


	alpha = A - B * g[1];
	v[0] = (F - B * v[1]) / alpha;

	// Calculate dU
	dU_new[i][0] = v[0];
	for (int j = 1; j < Ny; ++j) {
		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];
	}

}

void Solver::solve_middle_line_viscous(const int i) {

	Matrix Eti(4, Vector(4)), Ebi(4, Vector(4));

	static Matrix A(4, Vector(4)); // Relates to U(j+1) 
	static Matrix B(4, Vector(4)); // Relates to U(j) 
	static Matrix C(4, Vector(4)); // Relates to U(j-1) 
	static Vector F(4); // Right hand side  

	// Intermediate matrices
	static Matrix alpha(4, Vector(4));
	static Matrix v(Ny, Vector(4));
	static Tensor g(Ny, Matrix(4, Vector(4)));

	// Grab boundary values
	Ebi = inviscid_boundary_2D_E(BoundaryType.bottom, U[i][0], grid.jNorms(i, 0));
	Eti = inviscid_boundary_2D_E(BoundaryType.top, U[i][Ny - 1], grid.jNorms(i, Ny));
	

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity(4)
		- (i_minus_inviscid_Jacobians[i][Ny - 1] + i_viscous_Jacobians[i][Ny - 1]) * grid.iArea(i, Ny - 1)															// Left
		+ (i_plus_inviscid_Jacobians[i + 1][Ny - 1] - i_viscous_Jacobians[i + 1][Ny - 1]) * grid.iArea(i + 1, Ny - 1)												// Right
		- (j_minus_inviscid_Jacobians[i][Ny - 1] + j_viscous_Jacobians[i][Ny - 1]) * grid.jArea(i, Ny - 1)															// Bottom
		+ (j_plus_inviscid_Jacobians[i][Ny] + Eti * j_minus_inviscid_Jacobians[i][Ny] + j_viscous_Jacobians[i][Ny]) * grid.jArea(i, Ny);		// Right 

	C = (j_plus_inviscid_Jacobians[i][Ny - 1] - j_viscous_Jacobians[i][Ny - 1]) * (-grid.jArea(i, Ny - 1)); 

	F = (i_Fluxes[i][Ny - 1] + (i_plus_inviscid_Jacobians[i][Ny - 1] - i_viscous_Jacobians[i][Ny - 1]) * dU_old[i - 1][Ny - 1]) * (grid.iArea(i, Ny - 1)) 
		+ (i_Fluxes[i + 1][Ny - 1] + (i_minus_inviscid_Jacobians[i + 1][Ny - 1] + i_viscous_Jacobians[i + 1][Ny - 1]) * dU_old[i + 1][Ny - 1]) * (-grid.iArea(i + 1, Ny - 1))
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1)) 
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny)); 

	alpha = A;
	v[Ny - 1] = F / alpha;
	g[Ny - 1] = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		B = (j_minus_inviscid_Jacobians[i][j + 1] + j_viscous_Jacobians[i][j + 1]) * (grid.jArea(i, j + 1));

		A = grid.Volume(i, j) / dt * identity(4)
			- (i_minus_inviscid_Jacobians[i][j] + i_viscous_Jacobians[i][j]) * grid.iArea(i, j)
			+ (i_plus_inviscid_Jacobians[i + 1][j] - i_viscous_Jacobians[i + 1][j]) * grid.iArea(i + 1, j)
			- (j_minus_inviscid_Jacobians[i][j] + j_viscous_Jacobians[i][j]) * grid.jArea(i, j)
			+ (j_plus_inviscid_Jacobians[i][j + 1] - j_viscous_Jacobians[i][j + 1]) * grid.jArea(i, j + 1); 

		C = (j_plus_inviscid_Jacobians[i][j] - j_viscous_Jacobians[i][j]) * (-grid.jArea(i, j)); 

		F = (i_Fluxes[i][j] + (i_plus_inviscid_Jacobians[i][j] - i_viscous_Jacobians[i][j]) * dU_old[i - 1][j]) * (grid.iArea(i, j))
			+ (i_Fluxes[i + 1][j] + (i_minus_inviscid_Jacobians[i + 1][j] + i_viscous_Jacobians[i + 1][j]) * dU_old[i + 1][j]) * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j))
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1));	



		alpha = A - B * g[j + 1];
		g[j] = C / alpha;
		v[j] = (F - B * v[j + 1]) / alpha;
	}

	//  Bottom boundary
	B = (j_minus_inviscid_Jacobians[i][1] + j_viscous_Jacobians[i][1]) * (grid.jArea(i, 1));

	A = grid.Volume(i, 0) / dt * identity(4)
		- (i_minus_inviscid_Jacobians[i][0] + i_viscous_Jacobians[i][0]) * grid.iArea(i, 0)
		+ (i_plus_inviscid_Jacobians[i + 1][0] - i_viscous_Jacobians[i + 1][0]) * grid.iArea(i + 1, 0)
		- (Ebi * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0] + j_viscous_Jacobians[i][0]) * grid.jArea(i, 0) 
		+ (j_plus_inviscid_Jacobians[i][1] - j_viscous_Jacobians[i][1]) * grid.jArea(i, 1); 

	F = (i_Fluxes[i][0] + (i_plus_inviscid_Jacobians[i][0] - i_viscous_Jacobians[i][0]) * dU_old[i - 1][0]) * (grid.iArea(i, 0))
		+ (i_Fluxes[i + 1][0] + (i_minus_inviscid_Jacobians[i + 1][0] + i_viscous_Jacobians[i + 1][0]) * dU_old[i + 1][0]) * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1));


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
	Matrix Eti(4, Vector(4)), Ebi(4, Vector(4)), Eri(4, Vector(4));

	static Matrix A(4, Vector(4)); // Relates to U(j+1) 
	static Matrix B(4, Vector(4)); // Relates to U(j) 
	static Matrix C(4, Vector(4)); // Relates to U(j-1) 
	static Vector F(4); // Right hand side  

	// Intermediate matrices
	static Matrix alpha(4, Vector(4));
	static Matrix v(Ny, Vector(4));
	static Tensor g(Ny, Matrix(4, Vector(4)));

	// Grab boundary values
	Ebi = inviscid_boundary_2D_E(BoundaryType.bottom, U[i][0], grid.jNorms(i, 0));
	Eti = inviscid_boundary_2D_E(BoundaryType.top, U[i][Ny - 1], grid.jNorms(i, Ny)); 


	Eri = inviscid_boundary_2D_E(BoundaryType.right, U[i][Ny - 1], grid.iNorms(i + 1, Ny - 1));  

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity(4)
		- (i_minus_inviscid_Jacobians[i][Ny - 1] + i_viscous_Jacobians[i][Ny - 1]) * grid.iArea(i, Ny - 1)																// Left
		+ (i_plus_inviscid_Jacobians[i + 1][Ny - 1] + Eri * i_minus_inviscid_Jacobians[i + 1][Ny - 1] + i_viscous_Jacobians[i + 1][Ny - 1]) * grid.iArea(i + 1, Ny - 1)	// Right
		- (j_minus_inviscid_Jacobians[i][Ny - 1] + j_viscous_Jacobians[i][Ny - 1]) * grid.jArea(i, Ny - 1)																// Bottom
		+ (j_plus_inviscid_Jacobians[i][Ny] + Eti * j_minus_inviscid_Jacobians[i][Ny] + j_viscous_Jacobians[i][Ny]) * grid.jArea(i, Ny);		 						// Top 

	C = (j_plus_inviscid_Jacobians[i][Ny - 1] - j_viscous_Jacobians[i][Ny - 1]) * (-grid.jArea(i, Ny - 1)); 

	F = (i_Fluxes[i][Ny - 1] + (i_plus_inviscid_Jacobians[i][Ny - 1] - i_viscous_Jacobians[i][Ny - 1]) * dU_old[i - 1][Ny - 1]) * (grid.iArea(i, Ny - 1))
		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1))
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny));

	alpha = A;
	v[Ny - 1] = F / alpha;
	g[Ny - 1] = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		Eri = inviscid_boundary_2D_E(BoundaryType.right, U[i][j], grid.iNorms(i + 1, j));	

		B = (j_minus_inviscid_Jacobians[i][j + 1] + j_viscous_Jacobians[i][j + 1]) * (grid.jArea(i, j + 1)); 

		A = grid.Volume(i, j) / dt * identity(4)
			- (i_minus_inviscid_Jacobians[i][j] + i_viscous_Jacobians[i][j]) * grid.iArea(i, j)															// Left
			+ (i_plus_inviscid_Jacobians[i + 1][j] + Eri * i_minus_inviscid_Jacobians[i + 1][j] + i_viscous_Jacobians[i + 1][j]) * grid.iArea(i + 1, j)	// Right
			- (j_minus_inviscid_Jacobians[i][j] + j_viscous_Jacobians[i][j]) * grid.jArea(i, j)															// Bottom
			+ (j_plus_inviscid_Jacobians[i][j + 1] - j_viscous_Jacobians[i][j + 1]) * grid.jArea(i, j + 1);												// Top

		C = (j_plus_inviscid_Jacobians[i][j] - j_viscous_Jacobians[i][j]) * (-grid.jArea(i, j)); 

		F = (i_Fluxes[i][j] + (i_plus_inviscid_Jacobians[i][j] - i_viscous_Jacobians[i][j]) * dU_old[i - 1][j]) * (grid.iArea(i, j))
			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j))
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1)); 

		alpha = A - B * g[j + 1];
		g[j] = C / alpha;
		v[j] = (F - B * v[j + 1]) / alpha;
	}
	//  Bottom boundary
	Eri = inviscid_boundary_2D_E(BoundaryType.right, U[i][0], grid.iNorms(i + 1, 0));   

	B = (j_minus_inviscid_Jacobians[i][1] + j_viscous_Jacobians[i][1]) * (grid.jArea(i, 1)); 

	A = grid.Volume(i, 0) / dt * identity(4)
		- (i_minus_inviscid_Jacobians[i][0] + i_viscous_Jacobians[i][0]) * grid.iArea(i, 0)															// Left
		+ (i_plus_inviscid_Jacobians[i + 1][0] + Eri * i_minus_inviscid_Jacobians[i + 1][0] + i_viscous_Jacobians[i + 1][0]) * grid.iArea(i + 1, 0)	// Right 
		- (Ebi * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0] + j_viscous_Jacobians[i][0]) * grid.jArea(i, 0)					// Bottom 
		+ (j_plus_inviscid_Jacobians[i][1] - j_viscous_Jacobians[i][1]) * grid.jArea(i, 1);															// Top 

	F = (i_Fluxes[i][0] + (i_plus_inviscid_Jacobians[i][0] - i_viscous_Jacobians[i][0]) * dU_old[i - 1][0]) * (grid.iArea(i, 0))  
		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1));


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

	static Vector A(4); // Relates to U(j+1) 
	static Vector B(4); // Relates to U(j) 
	static Vector C(4); // Relates to U(j-1) 
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
	cout << inner_residual << endl;  
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

Vector Solver::minmod(const Vector& Ui, const Vector& Uii) { 
	
	Vector result(4);  

	for (int i = 0; i < 4; ++i) {
		result[i] = (Ui[i] * Uii[i] > 0) ? copysign(min(fabs(Ui[i]), fabs(Uii[i])), Ui[i]) : 0.0; 
	}

	return result; 
}

void Solver::write_2d_csv(const string& filename) {

	double a, density, pressure, u_velocity, v_velocity, Mach, Temperature;

	Vector Primitives(4);

	ofstream file(filename);
	file << "density, u velocity, v velocity, pressure, Mach, Temperature, x_points, y_points, z_points" << endl;

	 
	for (int i = 2; i < Nx + 2; ++i) { 
		for (int j = 2; j < Ny + 2; ++j) {  

			Primitives = constoPrim(U[i][j]);
			density = Primitives[0];
			pressure = Primitives[3];
			a = sqrt(gamma * pressure / density);
			u_velocity = Primitives[1];
			v_velocity = Primitives[2];
			Mach = sqrt(u_velocity * u_velocity + v_velocity * v_velocity) / a;
			Temperature = pressure / (density * R);

			file << density << ", " << u_velocity << ", " << v_velocity << ", " << pressure << ", " << Mach << ", " << Temperature<< ", " << grid.Center(i - 2, j - 2).x << ", " << grid.Center(i - 2, j - 2).y << ", 0.0" << endl;

		}
	}

	file.close();
	cout << "\033[36m2D File saved successfully as \"" << filename << "\"\033[0m" << endl;  
}

void Solver::write_1d_csv(const string& filename) {
	double a;

	Vector density(Ny), u_velocity(Ny), v_velocity(Ny), pressure(Ny), Mach(Ny), Temperature(Ny);
	Vector Primitives(Ny);


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

	file << "density, u velocity, pressure, Mach, Temperature, y_points" << endl;

	for (int j = 0; j < Ny; ++j) { 
		file << density[j]/INLET.rho << ", " << u_velocity[j]/INLET.u << ", "  << pressure[j]/INLET.p 
			<< ", " << Mach[j]/INLET.M << ", " << Temperature[j]/INLET.T << ", "  << grid.Center(Nx / 2, j).y << endl;  
	}
	

	file.close();
	cout << "\033[36m1D File saved successfully as \"" << filename << "\"\033[0m" << endl;

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////// Chemical Equilibrium ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


Chemistry::Chemistry(inlet_conditions& INLET) : rho(rho), u(u), v(v), T(T), h(h), h_g(h_g), h_chem(h_chem), h_stag(h_stag), p(p), EM(EM), CP_0(5), H_0(5), S_0(5), mu_k0(5), U(5), mk(5), nk(5), Yk(5), R_k(5),
		MW(5), theta_v(3), V(4), h_f(5), theta_f(5), R_mix(R_mix), T_coeff(5, Vector(7)), Int_const(5, Vector(2)), hi_t_coeff(5, Vector(7)), hi_t_int_const(5, Vector(2)), lo_t_coeff(5, Vector(7)),
		lo_t_int_const(5, Vector(2)), middle_t_coeff(5, Vector(7)), middle_t_int_const(5, Vector(2)) {

		EM.N2 = 0.78;
		EM.O2 = 1 - EM.N2;
		EM.total = EM.N2 + EM.O2;

		h_stag = cp * INLET.T + 0.5 * INLET.u * INLET.u;                    // Flow's stagnation enthalpy    
		MW = { 28.016, 32.0, 30.008, 14.008, 16 };                          // Molecular Weights of N2, O2, NO, N, O respectively

		for (int i = 0; i < 5; ++i) {
			R_k[i] = Ru / MW[i];
		}

		theta_v = { 3395.0, 2239.0, 2817.0 };                               // Characteristic temperatures of vibration for N2, O2, NO
		theta_f = { 0.0, 0.0, 2.996120e+6, 3.362160e+7, 1.542000e+7 };      // Enthalpies of formation 
		Yk = { EM.N2, EM.O2, 0, 0, 0 };                                     // Initialize for with mass fractions Yk


		// Three sets of NASA Polynomials for certain temperature ranges (200 - 1000, 1000 - 6000, 6000 - 20000)

		lo_t_coeff = {
			// N2, O2, NO, N, O
			{
				{2.210371497e+4, -3.818461820e+2, 6.082738360, -8.530914410e-3, 1.384646189E-05, -9.625793620e-9, 2.519705809e-12},
				{-3.425563420e+4, 4.847000970e+2, 1.119010961, 4.293889240e-3, -6.836300520e-7, -2.023372700e-9, 1.039040018e-12},
				{-1.143916503e+4, 1.536467592e+2, 3.431468730, -2.668592368e-3, 8.481399120e-6, -7.685111050e-9, 2.386797655e-12},
				{ 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0},
				{-7.953611300e+3, 1.607177787e+2, 1.966226438, 1.013670310e-3, -1.110415423e-6, 6.517507500e-10, -1.584779251e-13}
			}
		};

		lo_t_int_const = {
			{
				{7.108460860e+2, -1.076003744e+1},
				{-3.391454870e+3, 1.849699470e+1},
				{9.098214410e+3, 6.728725490},
				{5.610463780e+4, 4.193905036},
				{2.840362437e+4, 8.404241820}
			}
		};

		middle_t_coeff = {
			// N2, O2, NO, N, O
			{
				   { 5.877124060e+5, -2.239249073e+3, 6.066949220, -6.139685500e-4, 1.491806679e-7, -1.923105485e-11, 1.061954386e-15 },
				   { -1.037939022e+6, 2.344830282e+3, 1.819732036, 1.267847582e-3, -2.188067988e-7, 2.053719572e-11, -8.193467050e-16 },
				   { 2.239018716e+5, -1.289651623e+3, 5.433936030, -3.656034900e-4, 9.880966450e-8, -1.416076856e-11, 9.380184620e-16 },
				   { 8.876501380e+4, -1.071231500e+2, 2.362188287, 2.916720081e-4, -1.729515100e-7, 4.012657880e-11, -2.677227571e-15 },
				   { 2.619020262e+5, -7.298722030e+2, 3.317177270, -4.281334360e-4, 1.036104594e-7, -9.438304330e-12, 2.725038297e-16 }
			}
		};

		middle_t_int_const = {
			// N2, O2, NO, N, O     
			{
				{1.283210415e+4, -1.586640027e+1},
				{-1.689010929e+4, 1.738716506e+1},
				{1.750317656e+4, -8.501669090},
				{ 5.697351330e+4, 4.865231506 },
				{3.392428060e+4, -6.679585350e-1}
			}
		};

		hi_t_coeff = {
			// N2, O2, NO, N, O
			{
				{ 8.310139160e+8, -6.420733540e+5, 2.020264635e+2, -3.065092046e-2, 2.486903333e-6, -9.705954110e-11, 1.437538881e-15 },
				{ 4.975294300e+8, -2.866106874e+5, 6.690352250e+1, -6.169959020e-3, 3.016396027e-7, -7.421416600e-12, 7.278175770e-17 },
				{ -9.575303540e+8, 5.912434480e+5, -1.384566826e+2, 1.694339403e-2, -1.007351096e-6, 2.912584076e-11, -3.295109350e-16 },
				{ 5.475181050e+8, -3.107574980e+5, 6.916782740e+1, -6.847988130e-3, 3.827572400e-7, -1.098367709e-11, 1.277986024e-16 },
				{ 1.779004264e+8, -1.082328257e+5, 2.810778365e+1, -2.975232262e-3, 1.854997534e-7, -5.796231540e-12, 7.191720164e-17 }
			}
		};

		hi_t_int_const = {
			// N2, O2, NO, N, O
			{
				{4.938707040e+6, -1.672099740e+3},
				{2.293554027e+6, -5.530621610e+2},
				{-4.677501240e+6, 1.242081216e+3},
				{ 2.550585618e+6, -5.848769753e+2 },
				{8.890942630e+5, -2.181728151e+2}
			}
		};

		// Set polynomials based on temperature
		if (T > 200.0 && T < 1000.0) {
			T_coeff = lo_t_coeff;
			Int_const = lo_t_int_const;
		}
		else if (T >= 1000.0 && T < 6000.0) {
			T_coeff = middle_t_coeff;
			Int_const = middle_t_int_const;
		}
		else {
			T_coeff = hi_t_coeff;
			Int_const = hi_t_int_const;
		}

	}

void Chemistry::compute_h0() { 
	for (int j = 0; j < 5; ++j) {
		H_0[j] = Ru * T * (-T_coeff[j][0] * 1 / (T * T) + T_coeff[j][1] * log(T) / T + T_coeff[j][2] + T_coeff[j][3] * T / 2
			+ T_coeff[j][4] * T * T / 3 + T_coeff[j][5] * T * T * T / 4 + T_coeff[j][6] * T * T * T * T / 5 + Int_const[j][0] / T);
	}
}

void Chemistry::compute_s0() {
	for (int j = 0; j < 5; ++j) {
		S_0[j] = Ru * (-T_coeff[j][0] * 1 / (2 * T * T) - T_coeff[j][1] / T + T_coeff[j][2] * log(T) + T_coeff[j][3] * T
			+ T_coeff[j][4] * T * T / 2 + T_coeff[j][5] * T * T * T / 3 + T_coeff[j][6] * T * T * T * T / 4 + Int_const[j][1]);
	}
}

void Chemistry::compute_mass_fractions() {
	double sum = 0.0;
	for (int i = 0; i < 5; ++i) {
		mk[i] = nk[i] * MW[i];
		sum += mk[i];
	}


	for (int i = 0; i < 5; ++i) {
		Yk[i] = mk[i] / sum;
	}
}

void Chemistry::compute_mu0() {
	compute_h0(); 
	compute_s0(); 

	for (int j = 0; j < 5; ++j) {
		mu_k0[j] = H_0[j] - T * S_0[j];
	}
}

double Chemistry::norm(Vector& v1, Vector& v2) {
	double result = 0.0;
	for (int i = 0; i < 7; ++i) {
		result += fabs(v1[i] - v2[i]) * fabs(v1[i] - v2[i]);
	}
	return sqrt(result);
}

// Solve for molar fractions using Gibbs free energy minimization
void Chemistry::compute_molar_fractions() {
	compute_mu0();

	Vector X_new(8), X_old(8), dx(8); // Empty vectors

	// Initialize guesses
	X_old = { EM.total / 5, EM.total / 5, EM.total / 5, EM.total / 5, EM.total / 5, 0.0, 0.0, 0.0 };
	X_new = X_old / 2;

	double residual = norm(X_new, X_old);

	// Start equilibrium iterations
	while (residual >= 1e-10) {

		// set pi to 0 each iterations
		X_old[5] = 0; X_old[6] = 0;

		// Newton-Raphson Jacobian
		Matrix J = {
				{
				{1, 0, 0, 0, 0, -2, 0, -1},
				{0, 1, 0, 0, 0, 0, -2, -1},
				{0, 0, 1, 0, 0, -1, -1, -1},
				{0, 0, 0, 1, 0, -1, 0, -1},
				{0, 0, 0, 0, 1, 0, -1, -1},
				{2 * X_old[0], 0, X_old[2], X_old[3], 0, 0, 0, 0},
				{0, 2 * X_old[1], X_old[2], 0, X_old[4], 0, 0, 0},
				{X_old[0], X_old[1], X_old[2], X_old[3], X_old[4], 0, 0, -EM.total}
				}
		};

		// Gibbs free energy vector
		Vector G(5);
		for (int i = 0; i < 5; ++i) {
			if (X_old[i] <= 0) X_old[i] = 1e-10;
			G[i] = mu_k0[i] + Ru * T * log(X_old[i]) - Ru * T * log(EM.total) + Ru + T * log(p / 101325);
		}

		// Solution vector
		Vector F(8);
		F[0] = -G[0] / (Ru * T);
		F[1] = -G[1] / (Ru * T);
		F[2] = -G[2] / (Ru * T);
		F[3] = -G[3] / (Ru * T);
		F[4] = -G[4] / (Ru * T);
		F[5] = EM.N2 - (2 * X_old[0] + X_old[2] + X_old[3]);
		F[6] = EM.O2 - (2 * X_old[1] + X_old[2] + X_old[4]);
		F[7] = EM.total - (X_old[0] + X_old[1] + X_old[2] + X_old[3] + X_old[4]);

		// Update vector of unknowns 
		dx = F / J;

		// Solve for new molar fractions
		for (int i = 0; i < 8; ++i) {
			X_new[i] = X_old[i] * exp(dx[i]);
		}

		// Calculate residual and restart
		residual = norm(X_new, X_old);
		X_old = X_new;
	}

	// Final compositions
	for (int i = 0; i < 5; ++i) {
		nk[i] = X_new[i];
	}

	compute_mass_fractions();
};

void Chemistry:: display_mass_fraction() {
	cout << "Mass Fractions: " << endl << endl;
	cout << "N2: " << Yk[0] << endl;
	cout << "O2: " << Yk[1] << endl;
	cout << "NO: " << Yk[2] << endl;
	cout << "N:  " << Yk[3] << endl;
	cout << "O:  " << Yk[4] << endl;
	cout << "Temperature: " << T << " K" << endl;
}

// Solve for chemical equilibrium temperature using energy conservation
void Chemistry::compute_temperature() {

	double h_new = 0, cp_new = 0;
	T = h_stag / 1000;

	while (fabs(h_new - h) >= 0.1) {

		if (T > 200.0 && T < 1000.0) {
			T_coeff = lo_t_coeff;
			Int_const = lo_t_int_const;
		}
		else if (T >= 1000.0 && T < 6000.0) {
			T_coeff = middle_t_coeff;
			Int_const = middle_t_int_const;
		}
		else {
			T_coeff = hi_t_coeff;
			Int_const = hi_t_int_const;
		}

		compute_molar_fractions();

		h_new = 0.0;

		for (int i = 0; i < 3; ++i) {
			h_new += Yk[i] * (7 / 2 * R_k[i] * T + R_k[i] * theta_v[i] / (exp(theta_v[i] / T) - 1) + theta_f[i]);
		}

		for (int i = 3; i < 5; ++i) {
			h_new += Yk[i] * (5 / 2 * R_k[i] * T + theta_f[i]);
		}

		cp_new = h_new / T;

		T = T - 0.1 * (h_new - h) / cp_new;
	}
}

Vector Chemistry::compute_equilibrium(Vector& U) {
	Vector V = constoPrim(U);  // Convert conservative to primitive variables
	T = computeTemperature(U);

	R_mix = (EM.N2 / MW[0] + EM.O2 / MW[1]) * Ru; // Initial Mixture gas constant

	rho = V[0];
	u = V[1];
	v = V[2];
	p = rho * R_mix * T;  // Initial pressure estimate

	h = h_stag - 0.5 * sqrt(u * u + v * v);  // Static enthalpy

	compute_temperature();

	double e = 0.0;

	for (int i = 0; i < 3; ++i) {
		e += Yk[i] * (2.5 * R_k[i] * T + R_k[i] * (theta_v[i] / (exp(theta_v[i] / T) - 1))); 
	}

	for (int i = 3; i < 5; ++i) {
		e += Yk[i] * (1.5 * R_k[i] * T); 
	}

	// Recalculate R_mix with updated mass fractions
	R_mix = 0.0;
	for (int i = 0; i < 5; ++i) {
		R_mix += Yk[i] * R_k[i];
	}

	rho = p / (R_mix * T);
	
	double E = rho * e + 0.5 * rho * (u * u + v * v);
	U = { rho, rho * u, rho * v, E }; 


	return U;
}
