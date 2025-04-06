
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <fstream> 
#include "LinearAlgebra.h"
//#include "2DFVSLibrary.h" 

#include "2DFVSWithChemistry.h"

using namespace std; 


int main() {

	inlet_conditions INLET;

	INLET.p = 277,						// Inlet Pressure (SET)
	INLET.T = 251,						// Inlet Temperature (SET)
	INLET.M = 20,							// Inlet Mach speed (SET)
	INLET.a = 317.633,	// Inlet Sound Speed
	INLET.u = INLET.M * INLET.a,			// Inlet u-velocity
	INLET.v = 0,							// Inlet v-velocity
	INLET.rho = INLET.p / (R * INLET.T);	// Inlet density


	double rho, u, v, p, M, a;
	//M = 0.3804;
	//rho = 5.9259 * INLET.rho;
	//p = 466.4948 * INLET.p;  
	//a = sqrt(gamma * p / rho); 
	//u = M * a;  
	//v = 0.0;  

	M = INLET.M; 
	rho = INLET.rho; 
	p = INLET.p; 
	a = sqrt(gamma * p / rho); 
	u = M * a;  
	v = 0.0; 

	double gamma1 = 1.4;
	 

	Vector V = { rho, u, v, p };  
	Vector U = primtoCons(V, gamma1);   

	Chemistry chem(INLET);  

	Vector initial_moles = { 0.78, 0.22, 0.0, 0.0, 0.0 }; 
	Vector thermo = { 1.4, 287.0, 343.0 }; // gamma, R_mix
	displayVector(U); 

	cout << endl; 


	chem.compute_equilibrium(U, initial_moles, thermo);   
	chem.display_mass_fraction();

	displayVector(U); 

    return 0;
}