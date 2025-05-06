
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <fstream> 
#include "LinearAlgebra.h"
#include "GridGenerator.h" 
#include "2DFVSLibrary.h" 

using namespace std; 

int main() {

    int Nx = 50, Ny = 100; 

    CylinderGrid grid(Nx, Ny, 1, 3, 4.5, 0.0001, pi / 2, 3 * pi / 2);   

    int i = Nx - 1, j = Ny - 1;


	cout << grid.iNorms(i, j).x << " " << grid.iNorms(i, j).y << endl;  


    return 0;
}