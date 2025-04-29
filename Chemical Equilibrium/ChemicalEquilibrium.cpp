
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <fstream> 
#include "LinearAlgebra.h"
#include "2DFVSLibrary.h" 

using namespace std;

struct molar_concentration {
	double N2, O2, NO, N, O, total;
};

class Chemistry {
private:

	double rho, e, gamma, T, p, R_mix, cv_mix;
	Vector CP_0, H_0, S_0, mu_k0, U, mk, Xk, Yk, Ck, R_k, MW, theta_v, V, h_f, theta_f;
	Matrix T_coeff, Int_const, hi_t_coeff, hi_t_int_const, lo_t_coeff, lo_t_int_const, middle_t_coeff, middle_t_int_const;

	molar_concentration initial_moles;

public:

	Chemistry() : rho(rho), e(e), T(T), p(p), R_mix(R_mix), cv_mix(cv_mix), gamma(gamma),
		CP_0(5), H_0(5), S_0(5), mu_k0(5), U(5), mk(5), Xk(5), Yk(5), R_k(5), MW(5), theta_v(3), V(4), h_f(5), theta_f(5), Ck(5),
		T_coeff(5, Vector(7)), Int_const(5, Vector(2)), hi_t_coeff(5, Vector(7)), hi_t_int_const(5, Vector(2)), lo_t_coeff(5, Vector(7)),
		lo_t_int_const(5, Vector(2)), middle_t_coeff(5, Vector(7)), middle_t_int_const(5, Vector(2)) {

		// As a reminder, the order of species is N2, O2, NO, N, O, (Yk) is the mass fraction and (Xk) is the molar fraction. 

		MW = { 28.016, 32.0, 30.008, 14.008, 16 };	// Molecular Weights of N2, O2, NO, N, O respectively  

		// Calculate species gas constants
		for (int i = 0; i < 5; ++i) {
			R_k[i] = Ru / MW[i];
		}

		theta_v = { 3395.0, 2239.0, 2817.0 };                               // Characteristic temperatures of vibration for N2, O2, NO
		theta_f = { 0.0, 0.0, 2.996120e+6, 3.362160e+7, 1.542000e+7 };      // Enthalpies of formation 

		// Three sets of NASA Polynomials for certain temperature ranges [lo = 200 - 1000], middle = [1000 - 6000], hi = [6000 - 20000])
		lo_t_coeff = {
			// N2, O2, NO, N, O
			{
				{2.210371497e+4, -3.818461820e+2, 6.082738360, -8.530914410e-3, 1.384646189E-05, -9.625793620e-9, 2.519705809e-12},
				{-3.425563420e+4, 4.847000970e+2, 1.119010961, 4.293889240e-3, -6.836300520e-7, -2.023372700e-9, 1.039040018e-12},
				{-1.143916503e+4, 1.536467592e+2, 3.431468730, -2.668592368e-3, 8.481399120e-6, -7.685111050e-9, 2.386797655e-12},
				{ 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0},
				{-7.953611300e+3, 1.607177787e+2, 1.966226438, 1.013670310e-3, -1.110415423e-6, 6.517507500e-10, -1.584779251e-13}
			}
		};

		lo_t_int_const = {
			{
				{7.108460860e+2, -1.076003744e+1},
				{-3.391454870e+3, 1.849699470e+1},
				{9.098214410e+3, 6.728725490},
				{5.610463780e+4, 4.193905036},
				{2.840362437e+4, 8.404241820}
			}
		};

		middle_t_coeff = {
			// N2, O2, NO, N, O
			{
				   { 5.877124060e+5, -2.239249073e+3, 6.066949220, -6.139685500e-4, 1.491806679e-7, -1.923105485e-11, 1.061954386e-15 },
				   { -1.037939022e+6, 2.344830282e+3, 1.819732036, 1.267847582e-3, -2.188067988e-7, 2.053719572e-11, -8.193467050e-16 },
				   { 2.239018716e+5, -1.289651623e+3, 5.433936030, -3.656034900e-4, 9.880966450e-8, -1.416076856e-11, 9.380184620e-16 },
				   { 8.876501380e+4, -1.071231500e+2, 2.362188287, 2.916720081e-4, -1.729515100e-7, 4.012657880e-11, -2.677227571e-15 },
				   { 2.619020262e+5, -7.298722030e+2, 3.317177270, -4.281334360e-4, 1.036104594e-7, -9.438304330e-12, 2.725038297e-16 }
			}
		};

		middle_t_int_const = {
			// N2, O2, NO, N, O     
			{
				{1.283210415e+4, -1.586640027e+1},
				{-1.689010929e+4, 1.738716506e+1},
				{1.750317656e+4, -8.501669090},
				{ 5.697351330e+4, 4.865231506 },
				{3.392428060e+4, -6.679585350e-1}
			}
		};

		hi_t_coeff = {
			// N2, O2, NO, N, O
			{
				{ 8.310139160e+8, -6.420733540e+5, 2.020264635e+2, -3.065092046e-2, 2.486903333e-6, -9.705954110e-11, 1.437538881e-15 },
				{ 4.975294300e+8, -2.866106874e+5, 6.690352250e+1, -6.169959020e-3, 3.016396027e-7, -7.421416600e-12, 7.278175770e-17 },
				{ -9.575303540e+8, 5.912434480e+5, -1.384566826e+2, 1.694339403e-2, -1.007351096e-6, 2.912584076e-11, -3.295109350e-16 },
				{ 5.475181050e+8, -3.107574980e+5, 6.916782740e+1, -6.847988130e-3, 3.827572400e-7, -1.098367709e-11, 1.277986024e-16 },
				{ 1.779004264e+8, -1.082328257e+5, 2.810778365e+1, -2.975232262e-3, 1.854997534e-7, -5.796231540e-12, 7.191720164e-17 }
			}
		};

		hi_t_int_const = {
			// N2, O2, NO, N, O
			{
				{4.938707040e+6, -1.672099740e+3},
				{2.293554027e+6, -5.530621610e+2},
				{-4.677501240e+6, 1.242081216e+3},
				{ 2.550585618e+6, -5.848769753e+2 },
				{8.890942630e+5, -2.181728151e+2}
			}
		};

		// Set polynomials based on temperature
		if (T > 200.0 && T < 1000.0) {
			T_coeff = lo_t_coeff;
			Int_const = lo_t_int_const;
		}
		else if (T >= 1000.0 && T < 6000.0) {
			T_coeff = middle_t_coeff;
			Int_const = middle_t_int_const;
		}
		else {
			T_coeff = hi_t_coeff;
			Int_const = hi_t_int_const;
		}

	}

	void display_molar_fractions() {
		cout << "Molar Fractions: " << endl << endl;
		cout << "N2: " << Xk[0] << endl;
		cout << "O2: " << Xk[1] << endl;
		cout << "NO: " << Xk[2] << endl;
		cout << "N:  " << Xk[3] << endl;
		cout << "O:  " << Xk[4] << endl;
	}

	void display_mass_fraction() {
		cout << "Mass Fractions: " << endl << endl;
		cout << "N2: " << Yk[0] << endl;
		cout << "O2: " << Yk[1] << endl;
		cout << "NO: " << Yk[2] << endl;
		cout << "N:  " << Yk[3] << endl;
		cout << "O:  " << Yk[4] << endl;
		cout << "Temperature: " << T << endl;
		cout << "Pressure: " << p << endl;
		cout << "Internal energy: " << e << endl;
		cout << "Density: " << rho << endl;
		cout << "R_mix: " << R_mix << endl;
	}

	void compute_mass_fractions() {

		double MW_mix = 0.0; // Mixture Molecular Weight 

		for (int i = 0; i < 5; ++i) {
			MW_mix += Xk[i] * MW[i];
		}

		for (int i = 0; i < 5; ++i) {
			Yk[i] = (Xk[i] * MW[i]) / MW_mix;
		}
	}

	void compute_h0() {
		for (int j = 0; j < 5; ++j) {
			H_0[j] = Ru * T * (-T_coeff[j][0] * 1 / (T * T) + T_coeff[j][1] * log(T) / T + T_coeff[j][2] + T_coeff[j][3] * T / 2
				+ T_coeff[j][4] * T * T / 3 + T_coeff[j][5] * T * T * T / 4 + T_coeff[j][6] * T * T * T * T / 5 + Int_const[j][0] / T);
		}
	}

	void compute_s0() {
		for (int j = 0; j < 5; ++j) {
			S_0[j] = Ru * (-T_coeff[j][0] * 1 / (2 * T * T) - T_coeff[j][1] / T + T_coeff[j][2] * log(T) + T_coeff[j][3] * T
				+ T_coeff[j][4] * T * T / 2 + T_coeff[j][5] * T * T * T / 3 + T_coeff[j][6] * T * T * T * T / 4 + Int_const[j][1]);
		}
	}

	void compute_mu0() {
		compute_h0();
		compute_s0();

		for (int j = 0; j < 5; ++j) {
			mu_k0[j] = H_0[j] - T * S_0[j];
		}
	}

	double norm(Vector& v1, Vector& v2) {
		double result = 0.0;
		for (int i = 0; i < 7; ++i) {
			result += fabs(v1[i] - v2[i]) * fabs(v1[i] - v2[i]);
		}
		return sqrt(result);
	}

	bool safe_compute_molar_fractions() {

		try {
			compute_molar_fractions(); // Try to solve chemical equilibrium
			return true;
		}
		catch (const exception& ex) {
			cerr << "[Warning] Chemical equilibrium failed: " << ex.what() << std::endl;

			// Fallback: assume frozen air (N2 and O2 only)
			// N2 = 78%, O2 = 21%, everything else 0%
			Xk = { 0.78, 0.22, 0.0, 0.0, 0.0 };
			return false;
		}
	}

	// Solve for molar fractions using Gibbs free energy minimization
	void compute_molar_fractions() {
		compute_mu0(); // Recalculate reference chemical potentials

		Vector X_new(8), X_old(8), dx(8);

		// Initialize guess: equal distribution among species
		X_old = { initial_moles.total / 5, initial_moles.total / 5, initial_moles.total / 5, initial_moles.total / 5, initial_moles.total / 5, 0.0, 0.0, 0.0 };
		X_new = X_old / 2; // Slightly disturbed initial guess

		double residual = norm(X_new, X_old);
		int iteration = 0;
		const int max_iterations = 100; // Cap iterations to prevent infinite loops

		while (residual >= 1e-6 && iteration < max_iterations) {
			iteration++;

			// Reset intermediate variables
			X_old[5] = 0; X_old[6] = 0;

			// Build Jacobian
			Matrix J = {
				{
				{1, 0, 0, 0, 0, -2, 0, -1},
				{0, 1, 0, 0, 0, 0, -2, -1},
				{0, 0, 1, 0, 0, -1, -1, -1},
				{0, 0, 0, 1, 0, -1, 0, -1},
				{0, 0, 0, 0, 1, 0, -1, -1},
				{2 * X_old[0], 0, X_old[2], X_old[3], 0, 0, 0, 0},
				{0, 2 * X_old[1], X_old[2], 0, X_old[4], 0, 0, 0},
				{X_old[0], X_old[1], X_old[2], X_old[3], X_old[4], 0, 0, -initial_moles.total}
				}
			};

			// Gibbs free energy vector
			Vector G(5);
			for (int i = 0; i < 5; ++i) {
				double Xi_safe = max(X_old[i], 1e-10); // Protect log from zero
				G[i] = mu_k0[i] + Ru * T * log(Xi_safe / initial_moles.total) + Ru * T * log(p / 101325);
			}

			// Build residual vector
			Vector F(8);
			F[0] = -G[0] / (Ru * T);
			F[1] = -G[1] / (Ru * T);
			F[2] = -G[2] / (Ru * T);
			F[3] = -G[3] / (Ru * T);
			F[4] = -G[4] / (Ru * T);
			F[5] = initial_moles.N2 - (2 * X_old[0] + X_old[2] + X_old[3]);
			F[6] = initial_moles.O2 - (2 * X_old[1] + X_old[2] + X_old[4]);
			F[7] = initial_moles.total - (X_old[0] + X_old[1] + X_old[2] + X_old[3] + X_old[4]);

			// Solve for Newton update
			dx = F / J;

			// Update molar fractions safely
			for (int i = 0; i < 8; ++i) {
				double dx_safe = min(max(dx[i], -50.0), 50.0); // Clamp dx
				X_new[i] = X_old[i] * exp(dx_safe);

				if (std::isnan(X_new[i]) || std::isinf(X_new[i])) {
					throw std::runtime_error("Nonphysical molar fraction computed (NaN or Inf)");
				}
			}

			residual = norm(X_new, X_old);
			X_old = X_new;
		}

		if (iteration >= max_iterations) {
			throw std::runtime_error("Chemical equilibrium solver did not converge");
		}

		// Store final molar fractions
		for (int i = 0; i < 5; ++i) {
			Xk[i] = X_new[i];
		}
	}

	// Solve for chemical equilibrium temperature using energy conservation
	void compute_temperature() {

		double e_new = 0, cv_new = 0;

		while (fabs(e_new - e) >= 1) {

			if (T > 200.0 && T < 1000.0) {
				T_coeff = lo_t_coeff;
				Int_const = lo_t_int_const;
			}
			else if (T >= 1000.0 && T < 6000.0) {
				T_coeff = middle_t_coeff;
				Int_const = middle_t_int_const;
			}
			else {
				T_coeff = hi_t_coeff;
				Int_const = hi_t_int_const;
			}

			bool equilibrium_success = safe_compute_molar_fractions();
			compute_mass_fractions();

			if (equilibrium_success == false) {
				R_mix = 287.0;
				cv_mix = 717;
				T = e / cv_mix;
				p = rho * R_mix * T;
				break;
			}
			else {
				e_new = 0.0;

				for (int i = 0; i < 3; ++i) {
					e_new += Yk[i] * (2.5 * R_k[i] * T + R_k[i] * theta_v[i] / (exp(theta_v[i] / T) - 1) + theta_f[i]);
				}

				for (int i = 3; i < 5; ++i) {
					e_new += Yk[i] * (1.5 * R_k[i] * T + theta_f[i]);
				}

				cv_new = e_new / T;

				T = T - 0.1 * (e_new - e) / cv_new;

				R_mix = Yk[0] * R_k[0] + Yk[1] * R_k[1] + Yk[2] * R_k[2] + Yk[3] * R_k[3] + Yk[4] * R_k[4]; // Recalculate R_mix with updated mass fractions 
				p = rho * R_mix * T;
			}


		}
		cv_mix = cv_new;
	}

	Vector compute_equilibrium(double Rho, double E, Vector& initial_mol) {

		R_mix = 0.0;
		initial_moles.total = 0.0;

		initial_moles.N2 = initial_mol[0];
		initial_moles.O2 = initial_mol[1];
		initial_moles.NO = initial_mol[2];
		initial_moles.N = initial_mol[3];
		initial_moles.O = initial_mol[4];


		for (int i = 0; i < 5; ++i) {
			initial_moles.total += initial_mol[i];
			R_mix += (initial_mol[i] / MW[i]) * Ru;
		}

		rho = Rho;
		e = E;
		T = e / 717;
		p = rho * R_mix * T;  // Initial pressure estimate	  

		compute_temperature();
		gamma = 1 + R_mix / cv_mix;
		Vector thermo = { rho, e, p, T, R_mix, cv_mix, gamma, 0.0, 0.0 }; // Mass fractions    

		//////  Compute derivatives ////

		// dp/dp | e
		e = E + 1;
		T = e / 717;
		p = rho * 287.0 * T;
		compute_temperature();
		double p1 = p;


		e = E - 1;
		T = e / 717;
		p = rho * 287.0 * T;
		compute_temperature();
		double p2 = p;

		double dpde = (p1 - p2) / 2; // dp/de 


		// dp/de| rho
		rho = Rho + 0.00001;
		T = E / 717;
		p = rho * 287.0 * T;
		compute_temperature();
		p1 = p;
			
		rho = Rho - 0.00001;
		T = E / 717;
		p = rho * 287.0 * T;
		compute_temperature();
		p2 = p;
		double dpdrho = (p1 - p2) / 0.00002; // dp/drho

		thermo[7] = dpdrho; // dp/drho | e 
		thermo[8] = dpde; // dp/de | rho  

		return thermo;
	}



};

int main() {

	int n = 200;

	double e;
	double rho = 0.2;
	Chemistry chem;
	Vector thermo;
	Vector initial_moles = { 0.78, 0.22, 0.0, 0.0, 0.0 };
	ofstream file("Look_up_table.csv");
	file << "rho, e, p, T, R, cv, gamma, dpdrho, dpde" << endl; // Header for the CSV file  


	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			rho = 1e-4 + (10.0 - 1e-4) / n * i;
			e = 3e5 + (2e7 - 3e5) / n * j;
			thermo = chem.compute_equilibrium(rho, e, initial_moles);
			file << rho << ", " << e << ", " << thermo[2] << ", " << thermo[3] << ", " << thermo[4] << ", " << thermo[5] << ", " << thermo[6] << ", " << thermo[7] << ", " << thermo[8] << endl; // Write to CSV file  
		}
		cout << i * n << endl;
	}

	return 0;
}