
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <fstream> 
#include "GridGenerator.h" 
#include "2DFVSLibrary.h" 
#include "LinearAlgebra.h"
#include <omp.h> 
#include <array> 

using namespace std;



#define TIME chrono::high_resolution_clock::now(); 
#define DURATION chrono::duration<double> duration;   
// Function to multiply a fixed-size 2D array by a 1D array

array<double, 3> operator*(const array<array<double, 3>, 3>& A, const array<double, 3>& x) {
    array<double, 3> result = { 0, 0, 0 }; 

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}

int main() {
    array<array<double, 3>, 3> A = { {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}} };
    array<double, 3> x = { 1, 2, 3 };
    array<double, 3> result;

    auto start = TIME;  
    for (int i = 0; i < 100000; ++i) {
       result = A * x; // Uses overloaded operator*
    }
    auto end = TIME;

    DURATION duration = end - start;
    cout << duration.count() << endl; 


    Matrix B = random(3, 3);
    Vector C = random(3); 
    Vector D(3);

    start = TIME;
    for (int i = 0; i < 100000; ++i) {
        D = B * C;
    }
    end = TIME;

    duration = end - start;
    cout << duration.count() << endl;

    return 0;
}