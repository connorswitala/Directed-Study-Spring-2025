
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <fstream> 
#include <chrono>
#include "LinearAlgebra.h"
#include "2DFVSLibrary.h"
#include "GridGenerator.h"

using namespace std;


int main() {

	int n = 4, Nx = 5, Ny = 5;
	RampGrid grid(100, 100, 1, 1, 1, 0.75, 10); 

	auto start1 = chrono::high_resolution_clock::now(); 
 
	Tensor U = randomTensor(n, Nx, Ny);
	Tesseract
	Vector VL{ 0.15, 1.00, 0.00, 0.3 };
	Vector VR = { 0.02, 0.3, 0.00, 0.15 };

	Vector UL = primtoCons(VL);
	Vector UR = primtoCons(VR);

    U[0][0] = A_Plus(UL, grid.RNorms(0, 0)) * UL + A_Minus(UR, grid.RNorms(0, 0)) * UR;

	displayVector(U[0][0]);

	auto end1 = chrono::high_resolution_clock::now();
	chrono::duration<double> duration1 = end1 - start1;

	cout << endl << duration1.count(); 
	return 0;
} 