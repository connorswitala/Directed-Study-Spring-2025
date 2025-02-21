
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <fstream> 
#include "LinearAlgebra.h"
#include "2DFVSLibrary.h"
#include "GridGenerator.h"
#include <omp.h> 


using namespace std; 

int main() {

	SquareGrid grid(Nx, Ny, 1, 1, 0.01);

    cout << grid.jNorms(Nx/2 - 1, 0).x << endl;       


	string A = "plot_2D_grid.csv"; 
	ofstream file(A);  


	file << "x_points, y_points, z_points" << endl; 

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			file << grid.Center(i, j).x << ", " << grid.Center(i, j).y << ", 0.0" << endl; 
		}
	}

	file.close();

    return 0;
}