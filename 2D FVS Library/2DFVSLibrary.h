#pragma once

#include "LinearAlgebra.h"
#include "GridGenerator.h" 
#include <cstdlib>
#include <fstream> 
#include <omp.h>
#include <chrono> 
#include <iomanip>

#define TIME chrono::high_resolution_clock::now(); 
#define DURATION chrono::duration<double> duration; 

using namespace std;

// Global constants for fluid dynamics
constexpr double gamma = 1.4;
constexpr double R = 287.0;
constexpr double Ru = 8314.0; 
constexpr double cv = R / (gamma - 1);
constexpr double cp = cv + R;
constexpr double Pr = 0.71;

Point operator*(int s, const Point& normals); 



inline double computeInternalEnergy(const Vector& U) {
	return U[3] / U[0] - 1.0 / (2.0 * U[0] * U[0]) * (U[1] * U[1] + U[2] * U[2]);
} 
inline double computePressure(const Vector& U) {
	double e = computeInternalEnergy(U); 
	return (gamma - 1.0) * U[0] * e;  
}
inline double computeTemperature(const Vector& U) {
	double e = computeInternalEnergy(U); 
	return e / cv;  
}


// This enum class is for setting boundary conditions types
enum class BoundaryCondition {
	IsothermalWall,
	AdiabaticWall, 
	Inlet,
	Outlet,
	Symmetry,
	Undefined
};

// This struct contains the boundary conditions types for each side of the grid (left, right, bottom, top) 
struct BoundaryConditions {

	BoundaryCondition left;
	BoundaryCondition right;
	BoundaryCondition bottom;
	BoundaryCondition top;

	BoundaryConditions(BoundaryCondition left, BoundaryCondition right, BoundaryCondition bottom, BoundaryCondition top) : left(left), right(right), bottom(bottom), top(top) {}

};

// This set boundary conditions based on text from the UI
inline BoundaryCondition getBoundaryCondition(const string& input) {
	if (input == "inlet") return BoundaryCondition::Inlet;
	if (input == "outlet") return BoundaryCondition::Outlet;
	if (input == "symmetry") return BoundaryCondition::Symmetry;
	if (input == "adiabatic") return BoundaryCondition::AdiabaticWall;
	if (input == "isothermal") return BoundaryCondition::IsothermalWall;
	return BoundaryCondition::Undefined;  
}

// This sets the inlet flow conditions from inputs in the UI
struct inlet_conditions {
	double rho, u, v, p, T, M, a;
};

// This struct contains the states for inviscid Jacobian computation
struct Inviscid_State {
	double rho, u, v, a, p, uprime;  
};

// This struct contains the states for viscous Jacobians computation
struct Viscous_State { 
	double rho, u, v, p, a, k, uprime, pp, h0, T; 
};

// This function computes the states for the inviscid Jacobians
inline Inviscid_State compute_inviscid_state(const Vector& U, double nx, double ny) {
	Inviscid_State S;

	S.rho = U[0];
	S.u = U[1] / U[0]; 
	S.v = U[2] / U[0];  
	S.p = computePressure(U);
	S.a = sqrt(gamma * S.p / U[0]); 
	S.uprime = S.u * nx + S.v * ny;

	return S; 
}

// This function computes the states for the viscous Jacobians
inline Viscous_State compute_viscous_state(const Vector& U, double nx, double ny) { 
	Viscous_State S;
	double e = computeInternalEnergy(U);

	S.rho = U[0];
	S.u = U[1] / S.rho;
	S.v = U[2] / S.rho;
	S.p = computePressure(U);
	S.T = S.p / (S.rho * R);
	S.a = sqrt(gamma * S.p / S.rho);
	S.k = 1.0 / (S.a * S.a);
	S.uprime = S.u * nx + S.v * ny; 
	S.pp = (gamma - 1.0) * e; 
	S.h0 = U[3] / S.rho;
	return S;
}


Vector primtoCons(const Vector& V);
Vector constoPrim(const Vector& U);

class Solver { 

private:

	string gridtype; 

	const int Nx, Ny, progress_update; 
	double CFL, Tw, dt, inner_residual, t_tot;  

	Vector V_inlet, U_inlet, Global_Residual, t, iteration; 

	Tensor U, dU_new, dU_old, left_flux, right_flux, bot_flux, top_flux;
	Tesseract top_plus_jac, top_minus_jac, bottom_plus_jac, bottom_minus_jac, left_plus_jac, left_minus_jac, right_plus_jac, right_minus_jac; 

	Grid& grid;    
	BoundaryConditions BoundaryType;  
	inlet_conditions INLET; 


public: 

	double outer_residual;

	Solver(const int Nx, const int Ny, const inlet_conditions& INLET, Grid& grid, BoundaryConditions BoundaryType, double CFL, double Tw, int& progress_update);  

	Matrix inviscid_boundary_2D_E(BoundaryCondition type, const Vector& U, const Point& normals);
	Vector inviscid_boundary_2D_U(BoundaryCondition type, const Vector& U, const Point& normals);

	void get_ghost_cells();  
	void solve_inviscid();  
	void solve_inviscid_timestep(); 
	void compute_dt(); 
	void compute_inviscid_jacobians();

	void solve_left_line_inviscid();
	void solve_middle_line_inviscid(const int i);
	void solve_right_line_inviscid(); 


	void compute_inner_residual();
	void compute_outer_residual(); 
	void write_2d_csv(const string& filename); 
	void write_1d_csv(const string& filename);
	void write_residual_csv(); 
	
	void time(void (Solver::* func)()) { 
		auto start = TIME;

		(this->*func)();

		auto end = TIME;
		DURATION duration = end - start;
		cout << "Time taken: " << duration.count() << endl;
	}

	Vector minmod(Vector& Ui, Vector& Uii);

	//Vector constoViscPrim(const Vector& U);
	//Matrix viscous_boundary_2D_E(BoundaryCondition type, const Vector& U, const Point& normals);
	//Vector viscous_boundary_2D_U(BoundaryCondition type, const Vector& U, const Point& normals);
	//void solve_viscous();
	//void solve_viscous_timestep();
	//void compute_viscous_jacobians();
	//void solve_left_line_viscous();
	//void solve_middle_line_viscous(const int i);
	//void solve_right_line_viscous();
	//void viscous_calculations();

};






 