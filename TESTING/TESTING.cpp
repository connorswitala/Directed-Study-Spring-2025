
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

    RampGrid grid(Nx, Ny, 10, 10, 10, 10, 15);

    cout << grid.iNorms(50, 25).y << endl;

    return 0;
}