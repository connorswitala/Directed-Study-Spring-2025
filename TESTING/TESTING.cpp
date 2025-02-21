
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

Matrix outerProduct(const Vector& a, const Vector& b) {
    Matrix result = {}; // Initialize with zeros
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result[i][j] = a[i] * b[j];
        }
    }
    return result;
}

Vector Fluxes(const Vector& Vi, const Vector& Vii) {

	Vector Ui = primtoCons(Vi); 
	Vector Uii = primtoCons(Vii); 

	double nx = 1.0, ny = 0.0;
	double R = 287;
	double cv = R / (gamma - 1); 

	double g = 5.72; 
	double dp = fabs(Vii[3] - Vi[3]) / min(Vii[3], Vi[3]); 
	double weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1)); 
	Vector V = weight * Vi + (1 - weight) * Vii; 
	Vector U = primtoCons(V);   

	double rho = V[0], u = V[1], v = V[2], p = V[3]; 

	double uprime = u * nx + v * ny;
	double a = sqrt(gamma * p / rho);  

	double pe = (gamma - 1);
	double pp = 0.5 * (gamma - 1) * (u * u + v * v); 
	double h0 = (U[3] + p ) / rho;

	double l_plus = 0.5 * (uprime + fabs(uprime));
	double l_car_plus = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) + 0.5 * (uprime - a + fabs(uprime - a))) - l_plus; 
	double l_tild_plus = 0.5 * (0.5 * (uprime + a + fabs(uprime + a)) - 0.5 * (uprime - a + fabs(uprime - a)));


	Vector V1_PLUS = { l_car_plus / (a * a),
		(u * l_car_plus + a * nx * l_tild_plus) / (a * a),
		(v * l_car_plus + a * ny * l_tild_plus) / (a * a), 
		(h0 * l_car_plus + a * uprime * l_tild_plus) / (a * a) }; 

	Vector V2_PLUS = { l_tild_plus / a,
		u * l_tild_plus / a + nx * l_car_plus,
		v* l_tild_plus / a + ny * l_car_plus,
		h0* l_tild_plus / a + uprime * l_car_plus };  

	Vector m = { pp, -u * pe, -v * pe, pe }; 
	Vector n = { -uprime, nx, ny, 0 };

	Matrix A_Plus = outerProduct(V1_PLUS, m) + outerProduct(V2_PLUS, n) + l_plus * identity();

	double l_minus = 0.5 * (uprime - fabs(uprime)); 
	double l_car_minus = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) + 0.5 * (uprime - a - fabs(uprime - a))) - l_minus; 
	double l_tild_minus = 0.5 * (0.5 * (uprime + a - fabs(uprime + a)) - 0.5 * (uprime - a - fabs(uprime - a)));


	Vector V1_MINUS = { l_car_minus / (a * a),
	(u * l_car_minus + a * nx * l_tild_minus) / (a * a),
	(v * l_car_minus + a * ny * l_tild_minus) / (a * a),
	(h0 * l_car_minus + a * uprime * l_tild_minus) / (a * a) };

	Vector V2_MINUS = { l_tild_minus / a,
	u * l_tild_minus / a + nx * l_car_minus,
	v * l_tild_minus / a + ny * l_car_minus,
	h0 * l_tild_minus / a + uprime * l_car_minus };

	Matrix A_Minus = outerProduct(V1_MINUS, m) + outerProduct(V2_MINUS, n) + l_minus * identity(); 

	return A_Plus * Ui + A_Minus * Uii; 
}



int main() {

	Vector Vi = { 0.266, 0.927, 0.0, 0.303 };
	Vector Vii = { 0.125, 0.0, 0.0, 0.1 }; 
	Vector F; 
	auto start = TIME; 
	F = Fluxes(Vi, Vii);	
	auto end = TIME;
	DURATION duration = end - start;
	cout << "Fluxes took: " << duration.count() << " seconds." << endl;
	displayVector(F);

    return 0;
}