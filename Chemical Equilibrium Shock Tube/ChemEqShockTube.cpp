#include "LinearAlgebra.h"
#include "2DFVSWithChemistry.h"
#include "GridGenerator.h"

using namespace std; 


int main() {

	int Nx = 100, Ny = 100;

	Vector perfect_gas_air_molar_concentraiton = { 0.78, 0.22, 0.0, 0.0, 0.0 };


	Tensor Molar_Fractions(Nx, Matrix(Ny, Vector(5, 0.0)));   
	Tensor Thermo(Nx, Matrix(Ny, Vector(2, 0.0)));	// gamma then R_mix; 

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			Molar_Fractions[i][j] = perfect_gas_air_molar_concentraiton;
			Thermo[i][j] = { 1.4, 287.0 }; // gamma, R_mix 
		}
	}

	inlet_conditions INLET; 


	






}