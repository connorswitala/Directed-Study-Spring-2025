
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

	CylinderGrid grid(50, 25, 0.1, 0.3, 0.45, 0.001, pi / 2, 3 * pi / 2);

    cout << grid.jArea(45, 24) << endl;

    return 0;
}