#include <iostream>
#include <vector>
#include <cmath>
#include "LinearAlgebra.h"
#include "2DFVSLibrary.h"
#include "GridGenerator.h"


struct element_moles {
    double O2, N2, total;
};

class Chemistry {

private:

    double T, p, h_g, h, h_stag, R_mix, rho, u, v, h_chem; 
    Matrix T_coeff, Int_const, hi_t_coeff, hi_t_int_const, lo_t_coeff, lo_t_int_const, middle_t_coeff, middle_t_int_const;

    Vector Yk, mu_k0, H_0, S_0, CP_0, MW, mk, nk, theta_v, theta_f, V, U, h_f, R_k;  
    element_moles EM; 

public:

    Chemistry(element_moles& EM, inlet_conditions& INLET)       
        :rho(rho), u(u), v(v), T(T), h(h), h_g(h_g), h_chem(h_chem), h_stag(h_stag), p(p), EM(EM), CP_0(5), H_0(5), S_0(5), mu_k0(5), U(5), mk(5), nk(5), Yk(5), R_k(5),
        MW(5), theta_v(3), V(4), h_f(5), theta_f(5), R_mix(R_mix), T_coeff(5, Vector(7)), Int_const(5, Vector(2)), hi_t_coeff(5, Vector(7)), hi_t_int_const(5, Vector(2)), lo_t_coeff(5, Vector(7)),
        lo_t_int_const(5, Vector(2)), middle_t_coeff(5, Vector(7)), middle_t_int_const(5, Vector(2)) {

        h_stag = cp * INLET.T + 0.5 * INLET.u * INLET.u;                    // Flow's stagnation enthalpy    
        MW = { 28.016, 32.0, 30.008, 14.008, 16 };                          // Molecular Weights of N2, O2, NO, N, O respectively

        for (int i = 0; i < 5; ++i) {
            R_k[i] = Ru / MW[i]; 
        }

        theta_v = { 3395.0, 2239.0, 2817.0 };                               // Characteristic temperatures of vibration for N2, O2, NO
        theta_f = { 0.0, 0.0, 2.996120e+6, 3.362160e+7, 1.542000e+7 };      // Enthalpies of formation 
        Yk = { EM.N2, EM.O2, 0, 0, 0 };                                     // Initialize for with mass fractions Yk

         
        // Three sets of NASA Polynomials for certain temperature ranges (200 - 1000, 1000 - 6000, 6000 - 20000)

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

        void compute_h0() {
            for (int j = 0; j < 5; ++j) {
                H_0[j] =  Ru * T * (-T_coeff[j][0] * 1 / (T * T) + T_coeff[j][1] * log(T) / T + T_coeff[j][2] + T_coeff[j][3] * T / 2
                    + T_coeff[j][4] * T * T / 3 + T_coeff[j][5] * T * T * T / 4 + T_coeff[j][6] * T * T * T * T / 5 + Int_const[j][0] / T); 
            }
        }

        void compute_s0() { 
            for (int j = 0; j < 5; ++j) {
                S_0[j] = Ru * (-T_coeff[j][0] * 1 / (2 * T * T) - T_coeff[j][1] / T + T_coeff[j][2] * log(T) + T_coeff[j][3] * T 
                    + T_coeff[j][4] * T * T / 2 + T_coeff[j][5] * T * T * T / 3 + T_coeff[j][6] * T * T * T * T / 4 + Int_const[j][1]); 
            }
        }

        void compute_mass_fractions() {
            double sum = 0.0; 
            for (int i = 0; i < 5; ++i) {
                mk[i] = nk[i] * MW[i];  
                sum += mk[i];                  
            }


            for (int i = 0; i < 5; ++i) {
                Yk[i] = mk[i] / sum;              
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

        // Solve for molar fractions using Gibbs free energy minimization
        void compute_molar_fractions() {
            compute_mu0(); 
             
            Vector X_new(8), X_old(8), dx(8); // Empty vectors

            // Initialize guesses
            X_old = { EM.total / 5, EM.total / 5, EM.total / 5, EM.total / 5, EM.total / 5, 0.0, 0.0, 0.0 }; 
            X_new = X_old / 2;
            
            double residual = norm(X_new, X_old);   

            // Start equilibrium iterations
            while (residual >= 1e-10) { 

                // set pi to 0 each iterations
                X_old[5] = 0; X_old[6] = 0; 
                
                // Newton-Raphson Jacobian
                Matrix J = {
                        {
                        {1, 0, 0, 0, 0, -2, 0, -1}, 
                        {0, 1, 0, 0, 0, 0, -2, -1}, 
                        {0, 0, 1, 0, 0, -1, -1, -1}, 
                        {0, 0, 0, 1, 0, -1, 0, -1}, 
                        {0, 0, 0, 0, 1, 0, -1, -1}, 
                        {2 * X_old[0], 0, X_old[2], X_old[3], 0, 0, 0, 0},  
                        {0, 2 * X_old[1], X_old[2], 0, X_old[4], 0, 0, 0}, 
                        {X_old[0], X_old[1], X_old[2], X_old[3], X_old[4], 0, 0, -EM.total}
                        }
                };
                
                // Gibbs free energy vector
                Vector G(5);
                for (int i = 0; i < 5; ++i) {
                    if (X_old[i] <= 0) X_old[i] = 1e-10;   
                    G[i] = mu_k0[i] + Ru * T * log(X_old[i]) - Ru * T * log(EM.total) + Ru + T * log(p / 101325); 
                }

                // Solution vector
                Vector F(8);                
                F[0] = -G[0] / (Ru * T);   
                F[1] = -G[1] / (Ru * T);
                F[2] = -G[2] / (Ru * T);
                F[3] = -G[3] / (Ru * T);
                F[4] = -G[4] / (Ru * T);
                F[5] = EM.N2 - (2 * X_old[0] + X_old[2] + X_old[3]);   
                F[6] = EM.O2 - (2 * X_old[1] + X_old[2] + X_old[4]);   
                F[7] = EM.total - (X_old[0] + X_old[1] + X_old[2] + X_old[3] + X_old[4]); 

                // Update vector of unknowns 
                dx = F / J;  

                // Solve for new molar fractions
                for (int i = 0; i < 8; ++i) {
                    X_new[i] = X_old[i] * exp(dx[i]);
                }        

                // Calculate residual and restart
                residual = norm(X_new, X_old);   
                X_old = X_new; 
            }    

            // Final compositions
            for (int i = 0; i < 5; ++i) {
                nk[i] = X_new[i];     
            }

            compute_mass_fractions();         
        }; 

        void display_mass_fraction() {
            cout << "Mass Fractions: " << endl << endl;  
            cout << "N2: " << Yk[0] << endl;
            cout << "O2: " << Yk[1] << endl;
            cout << "NO: " << Yk[2] << endl;
            cout << "N:  " << Yk[3] << endl;
            cout << "O:  " << Yk[4] << endl;
            cout << "Temperature: " << T << " K" << endl;         
        }

        // Solve for chemical equilibrium temperature using energy conservation
        void compute_temperature() {
            
            double h_new = 0, cp_new = 0;
            T = h_stag / 1000; 

            while (fabs(h_new - h) >= 0.1) {              

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

                R_mix = 0.0; 
                for (int i = 0; i < 5; ++i) { 
                    R_mix += Yk[i] * R_k[i]; 
                } 
              
                p = rho * R_mix * T; 

                compute_molar_fractions();             

                h_new = 0.0;

                for (int i = 0; i < 3; ++i) {
                    h_new += Yk[i] * (7 / 2 * R_k[i] * T + R_k[i] * theta_v[i] / (exp(theta_v[i] / T) - 1) + theta_f[i]);  
                }

                for (int i = 3; i < 5; ++i) {
                    h_new += Yk[i] * (5 / 2 * R_k[i] * T + theta_f[i]);
                }

                cp_new = h_new / T;
                
                T = T - 0.1 * (h_new - h) / cp_new;      
            }              
        }


        Vector compute_equilibrium(Vector& U) {
            Vector V = constoPrim(U);  // Convert conservative to primitive variables
            T = computeTemperature(U);

            R_mix = (EM.N2 / MW[0] + EM.O2 / MW[1]) * Ru; // Initial Mixture gas constant

            rho = V[0];
            u = V[1];
            v = V[2];
            p = rho * R_mix * T;  // Initial pressure estimate

            h = h_stag - 0.5 * sqrt(u * u + v * v);  // Static enthalpy

            compute_temperature();

            // Recalculate R_mix with updated mass fractions
            R_mix = 0.0;
            for (int i = 0; i < 5; ++i) {
                R_mix += Yk[i] * R_k[i];
            }
             
            rho = p / (R_mix * T);
            double Z = p / (rho * R * T); 
            cout << "Z: " << Z << endl;


            return U;
        }


};

int main() {
    
    element_moles EM;
    EM.N2 = 0.78;  
    EM.O2 = 1 - EM.N2;    
    EM.total = EM.N2 + EM.O2;   
 
    inlet_conditions INLET;
    INLET.T = 251; 
    INLET.rho = 0.00385101; 
    INLET.p = INLET.rho * INLET.T * R;   
    INLET.a = sqrt(gamma * R * INLET.T);  
    INLET.M = 20;
    INLET.u = INLET.M * INLET.a; 
    INLET.v = 0; 


    inlet_conditions shock; 

    shock.M = 0.3804; 
    shock.p = 156560; 
    shock.T = 7000; 
    shock.rho = 0.05239;
    shock.a = sqrt(gamma * R * shock.T); 
    shock.u = shock.M * shock.a;  
    shock.v = 0; 

    //shock.M = 0.3804; 
    //shock.p = 466.4948 * INLET.p; 
    //shock.T = 7000;
    //shock.rho = 5.9259 * INLET.rho; 
    //shock.a = sqrt(gamma * R * shock.T); 
    //shock.u = shock.M * shock.a; 
    //shock.v = 0; 

    Vector V1 = { INLET.rho, INLET.u, INLET.v, INLET.p }; 
    Vector V = { shock.rho, shock.u, shock.v, shock.p }; 
    Vector U = primtoCons(V);  
    Vector U1 = primtoCons(V1);  

    Chemistry chem(EM, INLET);    


    Vector A = chem.compute_equilibrium(U);  
    chem.display_mass_fraction();


    return 0;
}
