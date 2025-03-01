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

constexpr double gamma = 1.4;
constexpr double R = 287.0;
constexpr double cv = R / (gamma - 1);
constexpr double cp = cv + R;
constexpr double Pr = 0.71;


// HERE //
constexpr int Nx = 100;
constexpr int Ny = 50;

using CellTensor = array < array < Vector, Ny>, Nx>;
using CellTesseract = array < array < Matrix, Ny>, Nx>;

using iFaceTensor = array < array < Vector, Ny>, Nx + 1>;
using iFaceTesseract = array < array <Matrix, Ny>, Nx + 1>;

using jFaceTensor = array < array < Vector, Ny + 1>, Nx>;
using jFaceTesseract = array < array <Matrix, Ny + 1>, Nx>;

typedef pair<Vector, Matrix> Duo;

enum class BoundaryCondition {
	IsothermalWall,
	AdiabaticWall, 
	Inlet,
	Outlet,
	Symmetry,
	Undefined
};

struct BoundaryConditions {

	BoundaryCondition left;
	BoundaryCondition right;
	BoundaryCondition bottom;
	BoundaryCondition top;

	BoundaryConditions(BoundaryCondition left, BoundaryCondition right, BoundaryCondition bottom, BoundaryCondition top) : left(left), right(right), bottom(bottom), top(top) {}

};

inline BoundaryCondition getBoundaryCondition(const string& input) {
	if (input == "inlet") return BoundaryCondition::Inlet;
	if (input == "outlet") return BoundaryCondition::Outlet;
	if (input == "symmetry") return BoundaryCondition::Symmetry;
	if (input == "adiabatic") return BoundaryCondition::AdiabaticWall;
	if (input == "isothermal") return BoundaryCondition::IsothermalWall;
	return BoundaryCondition::Undefined;  
}

struct inlet_conditions {
	double rho, u, v, p, T, M, a;
};

inline double computePressure(const Vector& U) {
	return (gamma - 1) * (U[3] - 0.5 * (U[1] * U[1] + U[2] * U[2])); 
}

inline double computeTemperature(const Vector& U) {
	return ((gamma - 1) * (U[3] - 0.5 * (U[1] * U[1] + U[2] * U[2]))) / (U[0] * R); 
}

struct Inviscid_State {
	double rho, u, v, p, a, k, uprime, pp, h0;
};

struct Viscous_State { 
	double rho, u, v, p, a, k, uprime, pp, h0, T; 
};

inline Inviscid_State compute_inviscid_state(const Vector& U, double nx, double ny) {
	Inviscid_State S;
	S.rho = U[0];
	S.u = U[1] / S.rho; 
	S.v = U[2] / S.rho;  
	S.p = computePressure(U); 
	S.a = sqrt(gamma * S.p / S.rho); 
	S.k = 1 / (S.a * S.a); 
	S.uprime = S.u * nx + S.v * ny;
	S.pp = 0.5 * (gamma - 1) * (S.u * S.u + S.v * S.v);
	S.h0 = (U[3] + S.p) / S.rho;
	return S; 
}

inline Viscous_State compute_viscous_state(const Vector& U, double nx, double ny) { 
	Viscous_State S;

	S.rho = U[0];
	S.u = U[1] / S.rho;
	S.v = U[2] / S.rho;
	S.p = computePressure(U);
	S.T = S.p / (S.rho * R);
	S.a = sqrt(gamma * S.p / S.rho);
	S.k = 1 / (S.a * S.a);
	S.uprime = S.u * nx + S.v * ny;
	S.pp = 0.5 * (gamma - 1) * (S.u * S.u + S.v * S.v);
	S.h0 = (U[3] + S.p) / S.rho;
	return S;
}


class Solver { 

private:

	CellTensor U, dU_new, dU_old;   

	iFaceTensor i_Fluxes; 
	iFaceTesseract i_plus_inviscid_Jacobians, i_minus_inviscid_Jacobians, i_viscous_Jacobians; 

	jFaceTensor j_Fluxes;
	jFaceTesseract j_plus_inviscid_Jacobians, j_minus_inviscid_Jacobians, j_viscous_Jacobians;

	Grid& grid;   
	double dt, inner_residual; 

	string gridtype; 

	BoundaryConditions BoundaryType; 
	double CFL, Tw; 
	Vector V_inlet, U_inlet;  
	inlet_conditions INLET; 
	int progress_update;  

public: 

	double outer_residual;

	Solver(const inlet_conditions& INLET, Grid& grid, BoundaryConditions BoundaryType, double CFL, double Tw, int& progress_update); 

	Vector primtoCons(const Vector& V);  
	Vector constoPrim(const Vector& U);
	Vector constoViscPrim(const Vector& U); 

	Matrix inviscid_boundary_2D_E(BoundaryCondition type, const Vector& U, const Point& normals);
	Vector inviscid_boundary_2D_U(BoundaryCondition type, const Vector& U, const Point& normals);

	Matrix viscous_boundary_2D_E(BoundaryCondition type, const Vector& U, const Point& normals);    
	Vector viscous_boundary_2D_U(BoundaryCondition type, const Vector& U, const Point& normals);   


	void solve_inviscid();  
	void solve_viscous();  
	void solve_inviscid_timestep(); 
	void solve_viscous_timestep(); 
	void compute_dt(); 
	void compute_inviscid_jacobians(); 
	void compute_viscous_jacobians(); 

	void solve_left_line_inviscid();
	void solve_middle_line_inviscid(const int i);
	void solve_right_line_inviscid(); 

	void solve_left_line_viscous();
	void solve_middle_line_viscous(const int i); 
	void solve_right_line_viscous(); 


	void compute_inner_residual();
	void compute_outer_residual(); 

	void write_2d_csv(const string& filename);
	void write_1d_csv(const string& filename);
	
	void time(void (Solver::* func)()) { 
		auto start = TIME;

		(this->*func)();

		auto end = TIME;
		DURATION duration = end - start;
		cout << "Time taken: " << duration.count() << endl;
	}


};






 