
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <fstream> 
#include "LinearAlgebra.h"
#include "GridGenerator.h" 
#include "2DFVSLibrary.h" 
#include <string>
#include <sstream>  // for std::stringstream

using namespace std; 

double find_closest(const std::vector<double>& subject, double target) {

    double low = subject[0]; 
    double high = subject[subject.size() - 1]; // Initialize last_value to the first element  
     
    if (target <= subject[low]) return low;   // target smaller than smallest 
    if (target >= subject[high]) return high; // target bigger than largest 

    while (low <= high) {
        int mid = (low + high) / 2;

        if (subject[mid] == target) {
            return mid;  // exact match
        }
        else if (subject[mid] < target) {
            low = mid + 1;
        }
        else {
            high = mid - 1;
        }
    }

    // Now low and high crossed over
    // Choose closer of array[low] and array[high] 
    if (low >= subject.size()) return high;
    if (high < 0) return low;

    if (std::abs(subject[low] - target) < std::abs(subject[high] - target)) {
        return low;
    }
    else {
        return high;
    }
}


int main() {

	Vector rho, e, R_mix, cv_mix, gamma, T, p, dpde, dpdrho;  

	ifstream file("C:/Users/Connor/source/repos/Directed Study Spring 2025/Chemical Equilibrium/Chemical_Equilibrium_Lookup_Table.csv"); 

	if (!file.is_open()) {
		std::cerr << "Could not open the file!\n";
		return 1;
	}

    string line;

    // Step 1: Read and ignore the header
    getline(file, line);

    // Step 2: Read the actual data
    while (getline(file, line)) {
        stringstream ss(line);
        string value;

        // Read each column value
        if (getline(ss, value, ',')) {
            rho.push_back(stod(value)); 
            cout << value << endl;
        }
        if (getline(ss, value, ',')) {
            e.push_back(stod(value));
        }
        if (getline(ss, value, ',')) {
            p.push_back(stod(value));
        }
        if (getline(ss, value, ',')) {
            T.push_back(stod(value)); 
        }
        if (getline(ss, value, ',')) {
            R_mix.push_back(stod(value)); 
        }
        if (getline(ss, value, ',')) {
            cv_mix.push_back(stod(value));  
        }
        if (getline(ss, value, ',')) {
            gamma.push_back(stod(value)); 
        }
        if (getline(ss, value, ',')) {
            dpdrho.push_back(stod(value)); 
        }
        if (getline(ss, value, ',')) {
            dpde.push_back(stod(value)); 
        } 
    }

    file.close();

    double target = 0.0254; 
	int index = find_closest(rho, target); // Find the index of the closest value in rho 
    cout << index;



}