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
	
	// IMPORTANT:  Due to how this code is set up, you can change Nx (number of cells in 'x' direction) and Ny (number of cells in 'y' direction) in the 
	// 2D FVS Library.h file. It should be right at the top. 

	auto start = TIME;

	double Wall_Temp = 300; // For Isothermal Wall boundary condition only (viscous solver). 
	double CFL = 1.0; // Pretty much only works at 1.0

	inlet_conditions INLET; 

	INLET.p = 100,						// Inlet Pressure (SET)
	INLET.T = 300,						// Inlet Temperature (SET)
	INLET.M = 15,							// Inlet Mach speed (SET)
	INLET.a = sqrt(perfGamma * perfR * INLET.T),	// Inlet Sound Speed
	INLET.u = INLET.M * INLET.a,			// Inlet u-velocity
	INLET.v = 0,							// Inlet v-velocity
	INLET.rho = INLET.p / (perfR * INLET.T);	// Inlet density 

	int progress_update = 1;  // This number prints a status update after the number of iterations declared here. 
	const int Nx = 100, Ny = 50; 

	BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Symmetry);      
	RampGrid grid(Nx, Ny, 10, 10, 10, 12, 25);    

	//BoundaryConditions BCs(BoundaryCondition::Outlet, BoundaryCondition::Outlet, BoundaryCondition::IsothermalWall, BoundaryCondition::Inlet);         
	//CylinderGrid grid(Nx, Ny, 0.1, 0.3, 0.45, 0.001, pi / 2, 3 * pi / 2); 

	//BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::IsothermalWall, BoundaryCondition::Symmetry);   
	//FlatPlateGrid grid(Nx, Ny, 1e-3, 1e-3, 5e-6);  

	//BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Symmetry);  
	//DoubleConeGrid grid(Nx, Ny, 1, 1, 1, 1, 25, 50, 2.5);  

	//BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Symmetry);   
	//MirroredGrid grid(Nx, Ny, 1, 1, 1, 2, 15, 30, 2.5);   
		

	// Sets up the solver class with inlet and boundary conditions, grid type, CFL, wall temperature, and the iteration number for progress     
	Solver solver(Nx, Ny, INLET, grid, BCs, CFL, Wall_Temp, progress_update); 

	// Calls solver type ( use -> solver.solve_viscous() if you want viscous solver)
	solver.solve_inviscid(); 

	printMemoryUsage(); // Prints maximum memory usage during the program.

	system("pause"); 
	return 0;
}



