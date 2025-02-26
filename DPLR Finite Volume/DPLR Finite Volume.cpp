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

	inlet_conditions INLET; 

	INLET.p = 10000.0,   
	INLET.T = 300.0,
	INLET.M = 2.5,  
	INLET.a = sqrt(gamma * R * INLET.T),   
	INLET.u = INLET.M * INLET.a, 
	INLET.v = 0,  
	INLET.rho = INLET.p / (R * INLET.T);   

	double Wall_Temp = 300; 

	const double CFL = 1.0; 

	

	BoundaryConditions BCs(BoundaryCondition::Outlet, BoundaryCondition::Outlet, BoundaryCondition::IsothermalWall, BoundaryCondition::Inlet);       
	CylinderGrid grid(Nx, Ny, 0.1, 0.3, 0.45, 0.00001, pi / 2, 3 * pi / 2); 

	//BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Symmetry);  
	//RampGrid grid(Nx, Ny, 10, 10, 10, 6, 15);

	//BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::IsothermalWall, BoundaryCondition::Symmetry);  
	//FlatPlateGrid grid(Nx, Ny, 1e-3, 1e-3, 5e-6);     

	const int progress_update = 1;  

	Solver solver(INLET, grid, BCs, CFL, Wall_Temp, progress_update);     

	solver.solve_viscous();       

	printMemoryUsage();

	return 0;
}