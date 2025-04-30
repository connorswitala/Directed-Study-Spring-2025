

#include <iostream>
#include <vector> 
using namespace std;

typedef vector<double> Vector;  

struct ThermoEntry {
    double rho, e, p, T, R, cv, gamma, dpdrho, dpde;
};


inline double computePressure(const Vector& U, double& gamma) {
    return (gamma - 1) * (U[3] - 0.5 * (U[1] * U[1] + U[2] * U[2]));
}

inline double computeSoundSpeed(const Vector& U, ThermoEntry& Thermo) {
    double p = computePressure(U, Thermo.gamma);
    return sqrt(p / (U[0] * U[0]) * Thermo.dpde + Thermo.dpdrho);
}

    int main() {


        double e = 500 * 717; 
        double rho = 0.3; 

        ThermoEntry test;
        test.cv = 717;
		test.gamma = 1.4; 
		test.R = 287;
        test.dpdrho = (test.gamma - 1) * e; 
		test.dpde = (test.gamma - 1) * rho; 
		test.p = (test.gamma - 1) * e * rho;   
		test.rho = rho;
		test.e = e;
		test.T = test.p / (test.rho * test.R);  
        
        Vector U = { rho, 0.0, 0.0, rho * (e) };
		double a = computeSoundSpeed(U, test);  

		cout << "Sound Speed: " << a << endl;

        return 0;
    }
