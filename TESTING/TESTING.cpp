
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

	Vector Ui, Uii, Vi, Vii, Up, Um, Flux;
	Matrix dudv, dvdu, M, Ap, Am;
	double nx, ny, l1, l2, l3, l4, pi, pii, dp, weight;
	Inviscid_State Si, Sii;

	Vi = { 0.15, 1.00, 0.00, 0.30 };
	Vii = { 0.02, 0.30, 0.00, 0.15 };

	// Left Face
	Ui = primtoCons(Vii);  
	Uii = primtoCons(Vi);  

	nx = -1;
	ny = 0;

	pi = computePressure(Ui);
	pii = computePressure(Uii);
	dp = fabs(pii - pi) / min(pii, pi);
	weight = 1;
	Up = weight * Ui + (1 - weight) * Uii;
	Um = (1 - weight) * Ui + weight * Uii;
	Si = compute_inviscid_state(Up, nx, ny);
	Sii = compute_inviscid_state(Um, nx, ny);

	l1 = 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a));
	l2 = 0.5 * (Si.uprime + fabs(Si.uprime));
	l3 = l2;
	l4 = 0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a));

	dudv = {
		{
			{1, 0, 0, 0},
			{Si.u, Si.rho, 0, 0},
			{Si.v, 0, Si.rho, 0},
			{0.5 * (Si.u * Si.u + Si.v * Si.v), Si.rho * Si.u, Si.rho * Si.v, cv / R}
		}
	};

	dvdu = {
	{
		{1, 0, 0, 0},
		{-Si.u / Si.rho, 1 / Si.rho, 0, 0},
		{-Si.v / Si.rho, 0, 1 / Si.rho, 0},
		{0.5 * R / cv * (Si.u * Si.u + Si.v * Si.v), -R / cv * Si.u, -R / cv * Si.v, R / cv}
	}
	};

	M = {
		{
			{l2, Si.rho * nx / (2 * Si.a) * (-l1 + l4), Si.rho * ny / (2 * Si.a) * (-l1 + l4), (l1 - 2 * l2 + l4) / (2 * Si.a * Si.a)},
			{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * Si.a * Si.rho) * (-l1 + l4) },
			{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * Si.a * Si.rho) * (-l1 + l4) },
			{0, 0.5 * Si.a * Si.rho * nx * (-l1 + l4), 0.5 * Si.a * Si.rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
		}
	};

	Ap = dudv * M * dvdu;

	l1 = 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a));
	l2 = 0.5 * (Sii.uprime - Sii.uprime);
	l3 = l2;
	l4 = 0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a));

	dudv = {
		{
			{1, 0, 0, 0},
			{Sii.u, Sii.rho, 0, 0},
			{Sii.v, 0, Sii.rho, 0},
			{0.5 * (Sii.u * Sii.u + Sii.v * Sii.v), Sii.rho * Sii.u, Sii.rho * Sii.v, cv / R}
		}
	};

	dvdu = {
		{
			{1, 0, 0, 0},
			{-Sii.u / Sii.rho, 1 / Sii.rho, 0, 0},
			{-Sii.v / Sii.rho, 0, 1 / Sii.rho, 0},
			{0.5 * R / cv * (Sii.u * Sii.u + Sii.v * Sii.v), -R / cv * Sii.u, -R / cv * Sii.v, R / cv}
		}
	};

	M = {
		{
			{l2, Sii.rho * nx / (2 * Sii.a) * (-l1 + l4), Sii.rho * ny / (2 * Sii.a) * (-l1 + l4), (l1 - 2 * l2 + l4) / (2 * Sii.a * Sii.a)},
			{0, 0.5 * (l1 * nx * nx + 2 * l3 * ny * ny + l4 * nx * nx), 0.5 * nx * ny * (l1 - 2 * l3 + l4), nx / (2 * Sii.a * Sii.rho) * (-l1 + l4) },
			{0, 0.5 * nx * ny * (l1 - 2 * l3 + l4), 0.5 * (l1 * ny * ny + 2 * l3 * nx * nx + l4 * ny * ny), ny / (2 * Sii.a * Sii.rho) * (-l1 + l4) },
			{0, 0.5 * Sii.a * Sii.rho * nx * (-l1 + l4), 0.5 * Sii.a * Sii.rho * ny * (-l1 + l4), 0.5 * (l1 + l4)}
		}
	};

	Am = dudv * M * dvdu;
	Flux = Ap * Ui + Am * Uii;

	displayVector(Flux); 

    return 0;
}