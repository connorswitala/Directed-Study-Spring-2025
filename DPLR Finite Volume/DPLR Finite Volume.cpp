#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <fstream> 
#include "GridGenerator.h" 
#include "2DFVSLibrary.h" 
#include "LinearAlgebra.h"
#include <iomanip>


using namespace std;


int main() {

	auto start = TIME;


	double p = 10000.0, T = 300.0, R = 287.0, M = 2.5, a = sqrt(gamma * R * T), u = M * a, v = 0, rho = p / (R * T); 
	const double CFL = 2.0; 
	Vector V_inlet = { rho, u, v, p };

	//BoundaryConditions BCs(BoundaryCondition::Outlet, BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Inlet);      
	//CylinderGrid grid(Nx, Ny, 0.1, 0.3, 0.45, 0.01, pi / 2, 3 * pi / 2); 


	BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Symmetry); 
	RampGrid grid(Nx, Ny, 10, 10, 10, 6, 15); 

	string gridtype; 
	if (dynamic_cast<RampGrid*>(&grid)) gridtype = "Ramp";
	else if (dynamic_cast<CylinderGrid*>(&grid)) gridtype = "Cylinder";
	else gridtype = "Unknown"; 

	cout << "Running DPLR for " << Nx << " by " << Ny << " " << gridtype << "..." << endl;

	Solver solver(V_inlet, grid, BCs, CFL);   

	solver.solve(); 

	string filename = to_string(Nx) + "x" + to_string(Ny) + "_" + gridtype + "_Solution.csv"; 
	solver.write2DCSV(filename);    

	printMemoryUsage();

	return 0;
}