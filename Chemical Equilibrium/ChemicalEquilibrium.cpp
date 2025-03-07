#include <iostream>
#include <vector>
#include <cmath>
#include "LinearAlgebra.h"

constexpr double R = 8.314;

struct element_moles {
    double O, N, total;
};

class Chemistry {

private:

    double T, P, h_g;  
    Matrix low_t_coeff, low_t_int_const, hi_t_coeff, hi_t_int_const;   

    Vector Yk, mu_k0, H_0, S_0, CP_0, MW, mk, nk, theta_vib;
    element_moles EM; 

public:

    Vector U; 

    Chemistry(double h_g, double T, double P, element_moles& EM)     
        : T(T), h_g(h_g), P(P), EM(EM), CP_0(5), H_0(5), S_0(5), mu_k0(5), U(5), mk(5), nk(5), Yk(5), MW(5), theta_vib(3) {         

        low_t_coeff = {
            // O, O2, N, N2, NO
            {
                { 2.619020262e+5, -7.298722030e+2, 3.317177270, -4.281334360e-4, 1.036104594e-7, -9.438304330e-12, 2.725038297e-16 },
                { -1.037939022e+6, 2.344830282e+3, 1.819732036, 1.267847582e-3, -2.188067988e-7, 2.053719572e-11, -8.193467050e-16 },
                { 8.876501380e+4, -1.071231500e+2, 2.362188287, 2.916720081e-4, -1.729515100e-7, 4.012657880e-11, -2.677227571e-15 },
                { 5.877124060e+5, -2.239249073e+3, 6.066949220, -6.139685500e-4, 1.491806679e-7, -1.923105485e-11, 1.061954386e-15 },
                { 2.239018716e+5, -1.289651623e+3, 5.433936030, -3.656034900e-4, 9.880966450e-8, -1.416076856e-11, 9.380184620e-16 }
            }
        };

        low_t_int_const = {
            // O, O2, N, N2, NO
            {
                {3.392428060e+4, -6.679585350e-1},
                {-1.689010929e+4, 1.738716506e+1},
                {5.697351330e+4, 4.865231506},
                {1.283210415e+4, -1.586640027e+1},
                {1.750317656e+4, -8.501669090}
            }
        };

        hi_t_coeff = {  
            // O, O2, N, N2, NO
            {
                { 1.779004264e+8, -1.082328257e+5, 2.810778365e+1, -2.975232262e-3, 1.854997534e-7, -5.796231540e-12, 7.191720164e-17 },
                { 4.975294300e+8, -2.866106874e+5, 6.690352250e+1, -6.169959020e-3, 3.016396027e-7, -7.421416600e-12, 7.278175770e-17 },
                { 5.475181050e+8, -3.107574980e+5, 6.916782740e+1, -6.847988130e-3, 3.827572400e-7, -1.098367709e-11, 1.277986024e-16 },
                { 8.310139160e+8, -6.420733540e+5, 2.020264635e+2, -3.065092046e-2, 2.486903333e-6, -9.705954110e-11, 1.437538881e-15 },
                { -9.575303540e+8, 5.912434480e+5, -1.384566826e+2, 1.694339403e-2, -1.007351096e-6, 2.912584076e-11, -3.295109350e-16 }
            }
        };

        hi_t_int_const = {
            // O, O2, N, N2, NO
            {
                {8.890942630e+5, -2.181728151e+2},
                {2.293554027e+6, -5.530621610e+2},
                {2.550585618e+6, -5.848769753e+2},
                {4.938707040e+6, -1.672099740e+3},
                {-4.677501240e+6, 1.242081216e+3}
            }
        };


        MW = { 15.999, 31.999, 14.0067, 28.0134, 44.013 }; 

        theta_vib = { 2239.0, 3352.0, 2690 }; 
    }

        void compute_cp0_low_t() { 
            for (int j = 0; j < 5; ++j) {
                
                CP_0[j] =  R * (low_t_coeff[j][0] * 1 / (T * T) + low_t_coeff[j][1] / T + low_t_coeff[j][2] + low_t_coeff[j][3] * T + low_t_coeff[j][4] * T * T 
                    + low_t_coeff[j][5] * T * T * T + low_t_coeff[j][6] * T * T * T * T);  
            }
        }

        void compute_h0_low_t() {
            for (int j = 0; j < 5; ++j) {
                H_0[j] =  R * T * (-low_t_coeff[j][0] * 1 / (T * T) + low_t_coeff[j][1] * log(T) / T + low_t_coeff[j][2] + low_t_coeff[j][3] * T / 2
                    + low_t_coeff[j][4] * T * T / 3 + low_t_coeff[j][5] * T * T * T / 4 + low_t_coeff[j][6] * T * T * T * T / 5 + low_t_int_const[j][0] / T); 
            }
        }

        void compute_s0_low_t() { 
            for (int j = 0; j < 5; ++j) {
                S_0[j] = R * (-low_t_coeff[j][0] * 1 / (2 * T * T) - low_t_coeff[j][1] / T + low_t_coeff[j][2] * log(T) + low_t_coeff[j][3] * T
                    + low_t_coeff[j][4] * T * T / 2 + low_t_coeff[j][5] * T * T * T / 3 + low_t_coeff[j][6] * T * T * T * T / 4 + low_t_int_const[j][1]);
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

        void compute_mu0_low_t() {  
            compute_h0_low_t();
            compute_s0_low_t();   
          
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

        void compute_molar_fractions_low_t() {
            compute_mu0_low_t();
             
            Vector X_new(8), X_old(8), dx(8); 
            X_old = { 0.5, 0.1, 0.2, 1, 0.2, 0.0, 0.0, 0.0 };   
            X_new = 10 * X_old;  
            
            double residual = norm(X_new, X_old); 
            
              
            while (residual >= 1e-6) { 

                X_old[5] = 0; X_old[6] = 0; 
                
                Matrix J = {
                        {
                        {1, 0, 0, 0, 0, -1, 0, -1}, 
                        {0, 1, 0, 0, 0, -2, 0, -1}, 
                        {0, 0, 1, 0, 0, 0, -1, -1}, 
                        {0, 0, 0, 1, 0, 0, -2, -1}, 
                        {0, 0, 0, 0, 1, -1, -1, -1}, 
                        {X_old[0], 2 * X_old[1], 0, 0, X_old[4], 0, 0, 0}, 
                        {0, 0, X_old[2], 2 * X_old[3], X_old[4], 0, 0, 0},
                        {X_old[0], X_old[1], X_old[2], X_old[3], X_old[4], 0, 0, -EM.total}
                        }
                };
           
                Vector G(5);

                for (int i = 0; i < 5; ++i) {
                    if (X_old[i] <= 0) X_old[i] = 1e-10;   
                    G[i] = mu_k0[i] + R * T * log(X_old[i]) - R * T * log(EM.total) + R + T * log(P / 101325); 
                }

                Vector F(8);
                
                F[0] = -G[0] / (R * T);   
                F[1] = -G[1] / (R * T);
                F[2] = -G[2] / (R * T);
                F[3] = -G[3] / (R * T);
                F[4] = -G[4] / (R * T);
                F[5] = EM.O - (X_old[0] + 2 * X_old[1] + X_old[4]); 
                F[6] = EM.N - (X_old[2] + 2 * X_old[3] + X_old[4]);  
                F[7] = EM.total - (X_old[0] + X_old[1] + X_old[2] + X_old[3] + X_old[4]); 
                
                dx = F / J;  
                for (int i = 0; i < 8; ++i) {
                    X_new[i] = X_old[i] * exp(dx[i]);
                }             
                residual = norm(X_new, X_old);   

                X_old = X_new; 
            }    

            for (int i = 0; i < 5; ++i) {
                nk[i] = (X_new[i]);     
            }

            compute_mass_fractions();  
       
        }; 

        void display_mass_fraction() {
            cout << "Mass Fractions: " << endl << endl;  
            cout << "O:   " << Yk[0] << endl;
            cout << "O2:  " << Yk[1] << endl;
            cout << "N:   " << Yk[2] << endl;
            cout << "N2:  " << Yk[3] << endl;
            cout << "NO:  " << Yk[4] << endl;
            cout << "Temperature: " << T << " K" << endl; 
        
        }

        void compute_cp0_hi_t() {
            for (int j = 0; j < 5; ++j) {

                CP_0[j] = R * (hi_t_coeff[j][0] * 1 / (T * T) + hi_t_coeff[j][1] / T + hi_t_coeff[j][2] + hi_t_coeff[j][3] * T + hi_t_coeff[j][4] * T * T
                    + hi_t_coeff[j][5] * T * T * T + hi_t_coeff[j][6] * T * T * T * T);
            }
        }

        void compute_h0_hi_t() {
            for (int j = 0; j < 5; ++j) {
                H_0[j] = R * T * (-hi_t_coeff[j][0] * 1 / (T * T) + hi_t_coeff[j][1] * log(T) / T + hi_t_coeff[j][2] + hi_t_coeff[j][3] * T / 2
                    + hi_t_coeff[j][4] * T * T / 3 + hi_t_coeff[j][5] * T * T * T / 4 + hi_t_coeff[j][6] * T * T * T * T / 5 + hi_t_int_const[j][0] / T);
            }
        }

        void compute_s0_hi_t() {
            for (int j = 0; j < 5; ++j) {
                S_0[j] = R * (-hi_t_coeff[j][0] * 1 / (2 * T * T) - hi_t_coeff[j][1] / T + hi_t_coeff[j][2] * log(T) + hi_t_coeff[j][3] * T
                    + hi_t_coeff[j][4] * T * T / 2 + hi_t_coeff[j][5] * T * T * T / 3 + hi_t_coeff[j][6] * T * T * T * T / 4 + hi_t_int_const[j][1]);
            }
        }

        void compute_mu0_hi_t() {
            compute_h0_hi_t();
            compute_s0_hi_t();

            for (int j = 0; j < 5; ++j) {
                mu_k0[j] = H_0[j] - T * S_0[j];
            }
        }

        void compute_molar_fractions_hi_t() {
            compute_mu0_hi_t(); 

            Vector X_new(8), X_old(8), dx(8);
            X_old = { EM.total/4, EM.total/1000, EM.total/4, EM.total/2, EM.total/1000, 0.0, 0.0, 0.0 }; 
            X_new = 10 * X_old;

            double residual = norm(X_new, X_old);


            while (residual >= 1e-8) {

                X_old[5] = 0; X_old[6] = 0;

                Matrix J = {
                        {
                        {1, 0, 0, 0, 0, -1, 0, -1},
                        {0, 1, 0, 0, 0, -2, 0, -1},
                        {0, 0, 1, 0, 0, 0, -1, -1},
                        {0, 0, 0, 1, 0, 0, -2, -1},
                        {0, 0, 0, 0, 1, -1, -1, -1},
                        {X_old[0], 2 * X_old[1], 0, 0, X_old[4], 0, 0, 0},
                        {0, 0, X_old[2], 2 * X_old[3], X_old[4], 0, 0, 0},
                        {X_old[0], X_old[1], X_old[2], X_old[3], X_old[4], 0, 0, -EM.total}
                        }
                };

                Vector G(5);

                for (int i = 0; i < 5; ++i) {
                    if (X_old[i] <= 0) X_old[i] = 1e-10;
                    G[i] = mu_k0[i] + R * T * log(X_old[i]) - R * T * log(EM.total) + R + T * log(P / 101325);  
                }

                Vector F(8);

                F[0] = -G[0] / (R * T);
                F[1] = -G[1] / (R * T);
                F[2] = -G[2] / (R * T);
                F[3] = -G[3] / (R * T);
                F[4] = -G[4] / (R * T);
                F[5] = EM.O - (X_old[0] + 2 * X_old[1] + X_old[4]);
                F[6] = EM.N - (X_old[2] + 2 * X_old[3] + X_old[4]);
                F[7] = EM.total - (X_old[0] + X_old[1] + X_old[2] + X_old[3] + X_old[4]);

                dx = F / J;
                for (int i = 0; i < 8; ++i) {
                    X_new[i] = X_old[i] * exp(dx[i]);
                }
                residual = norm(X_new, X_old);

                X_old = X_new;
            }

            for (int i = 0; i < 5; ++i) {
                nk[i] = (X_new[i]);
            }

            compute_mass_fractions();

        }; 

        void compute_temperature_hi_t() {
            
            double T_new = 0.9 * T, h = 0.999 * h_g, cp;   

            while (fabs(h - h_g) >= 1e-6) { 
                compute_molar_fractions_hi_t(); 
                compute_cp0_hi_t(); 

                h = 0.0, cp = 0.0; 

                for (int i = 0; i < 5; ++i) { 
                    h += Yk[i] * H_0[i]; 
                    cp += Yk[i] * CP_0[i];  
                }

                cout << h << endl; 
                
                
                T_new = T - (h - h_g) / cp;                  
                T = T_new;  
                cout << T << endl;
            }
        }

};

int main() {
    
    element_moles EM;
    EM.N = 0.767083; 
    EM.O = 1 - EM.N;  
    EM.total = EM.N + EM.O;   
 

    double h_g = 1000000;
    double P = 101325;
    double T_CFD = 800;  

    Chemistry chem(h_g, T_CFD, P, EM);   


    chem.compute_temperature_hi_t();
    chem.display_mass_fraction(); 


    return 0;
}
