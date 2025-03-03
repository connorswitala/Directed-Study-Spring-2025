
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <fstream> 
#include "LinearAlgebra.h"



using namespace std; 

int main() {

    Vector A = random(4), B = random(4);
    Matrix C = random(4, 4), D = random(4, 4); 


    Vector E = A/C;
    Matrix F = C/D;


    displayVector(E);
    displayMatrix(F); 

    return 0;
}