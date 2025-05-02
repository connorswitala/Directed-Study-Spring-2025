

#include <iostream>
#include <vector> 
#include <cstdlib>
#include <fstream> 
#include <iomanip>
#include <sstream>
#include <algorithm>

using namespace std;

typedef vector<double> Vector;  
typedef vector<vector<double>> Matrix;  

struct ThermoEntry {
    double rho, e, p, T, R, cv, gamma, dpdrho, dpde;
};


vector<vector<ThermoEntry>> load_csv(const string& filename) {
	ifstream file(filename);
	string line;
	vector<vector<ThermoEntry>> data(200, vector<ThermoEntry>(200));

	getline(file, line); // skip header

	int i = 0, j = 0;
	while (getline(file, line)) {
		stringstream ss(line);
		ThermoEntry entry;
		string temp;
		getline(ss, temp, ','); entry.rho = stod(temp);
		getline(ss, temp, ','); entry.e = stod(temp);
		getline(ss, temp, ','); entry.p = stod(temp);
		getline(ss, temp, ','); entry.T = stod(temp);
		getline(ss, temp, ','); entry.R = stod(temp);
		getline(ss, temp, ','); entry.cv = stod(temp);
		getline(ss, temp, ','); entry.gamma = stod(temp);
		getline(ss, temp, ','); entry.dpdrho = stod(temp);
		getline(ss, temp, ','); entry.dpde = stod(temp);

		data[i][j] = entry;

		// Advance j (energy index), then i (rho index)
		j++;
		if (j == 200) {
			j = 0;
			i++;
			if (i == 200) break; // just in case file has more lines
		}
	}
	return data;
}
int find_index(double target, double min, double max, int n) {
	double delta = (max - min) / (n);
	int idx = static_cast<int>((target - min) / delta);
	if (idx < 0) idx = 0;
	if (idx > n - 2) idx = n - 2;
	return idx;
}
ThermoEntry bilinear_interpolate(const std::vector<vector<ThermoEntry>>& table, double rho, double e) {

	int n_rho = 200, n_e = 200;
	double rho_min = 1e-4, rho_max = 10.0, e_min = 3e5, e_max = 2e7;

	int i = find_index(rho, rho_min, rho_max, n_rho);
	int j = find_index(e, e_min, e_max, n_e);

	double drho = (rho_max - rho_min) / (n_rho);
	double de = (e_max - e_min) / (n_e);

	double rho_i = rho_min + i * drho;
	double e_j = e_min + j * de;
	double t = (rho - rho_i) / drho;
	double u = (e - e_j) / de;

	auto get = [&](int ii, int jj) -> const ThermoEntry& {
		return table[ii][jj];
		};

	ThermoEntry Q11 = get(i, j);
	ThermoEntry Q21 = get(i + 1, j);
	ThermoEntry Q12 = get(i, j + 1);
	ThermoEntry Q22 = get(i + 1, j + 1);

	ThermoEntry result;

	auto lerp2D = [&](double a, double b, double c, double d) {
		return (1 - t) * (1 - u) * a + t * (1 - u) * b + (1 - t) * u * c + t * u * d;
		};

	result.rho = rho;
	result.e = e;
	result.T = lerp2D(Q11.T, Q21.T, Q12.T, Q22.T);
	result.R = lerp2D(Q11.R, Q21.R, Q12.R, Q22.R);
	result.cv = lerp2D(Q11.cv, Q21.cv, Q12.cv, Q22.cv);
	result.gamma = lerp2D(Q11.gamma, Q21.gamma, Q12.gamma, Q22.gamma);
	result.p = lerp2D(Q11.p, Q21.p, Q12.p, Q22.p);
	result.dpdrho = lerp2D(Q11.dpdrho, Q21.dpdrho, Q12.dpdrho, Q22.dpdrho);
	result.dpde = lerp2D(Q11.dpde, Q21.dpde, Q12.dpde, Q22.dpde);

	return result;
}


    int main() {

		string chemistry_table = "C:/Users/Connor/source/repos/Directed Study Spring 2025/Chemical Equilibrium/Chemical_Equilibrium_Lookup_Table.csv";  // path to your CSV file
		vector<vector<ThermoEntry>>cell_thermo = load_csv(chemistry_table); 

		double p = 10000;
		double T = 300;
		double R = 287.058;
		double M = 2.5;
		double a = sqrt(1.4 * R * T);
		double u = M * a;
		double v = 0;
		double rho = p / (R * T);

		double e = (1.4 - 1) * p / rho; // Internal energy
		cout << e << endl;


        return 0;
    }
