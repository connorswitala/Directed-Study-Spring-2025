#include <iostream>
#include <vector>
#include <cmath>
#include "LinearAlgebra.h"
#include "2DFVSLibrary.h"
#include "GridGenerator.h"



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

    displayVector(U); 
    
    Vector A = chem.compute_equilibrium(U);  

    displayVector(A); 



    return 0;
}
