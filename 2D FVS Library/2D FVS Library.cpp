// 2D FVS Library.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "framework.h"
#include "2DFVSLibrary.h"


Vector primtoCons(const Vector& V, double gamma) { 

	Vector U(4);
	U[0] = V[0];
	U[1] = V[0] * V[1];
	U[2] = V[0] * V[2];
	U[3] = V[3] / (gamma - 1) + 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2]);

	return U;
}
Vector constoPrim(const Vector& U, double gamma) { 

	static Vector V(4);
	V[0] = U[0];
	V[1] = U[1] / U[0];
	V[2] = U[2] / U[0];
	V[3] = (U[3] - 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2])) * (gamma - 1);

	return V;
}



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

	int n_rho = 400, n_e = 400;
	double rho_min = 1e-4, rho_max = 10.0, e_min = 3e5, e_max = 2e7;

	int i = find_index(rho, rho_min, rho_max, n_rho);
	int j = find_index(e, e_min, e_max, n_e);

	double drho = (rho_max - rho_min) / (n_rho - 1); 
	double de = (e_max - e_min) / (n_e - 1); 

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



Solver::Solver(const int Nx, const int Ny, const inlet_conditions& INLET, Grid& grid, BoundaryConditions BoundaryType, double CFL, double Tw, int& progress_update) : Nx(Nx), Ny(Ny), 
Tw(Tw), INLET(INLET), Global_Residual(Global_Residual), t(t), t_tot(t_tot), iteration(iteration), V_inlet(V_inlet), U_inlet(U_inlet), grid(grid), BoundaryType(BoundaryType), U(U), dU_new(dU_new), dU_old(dU_old), i_Fluxes(i_Fluxes), j_Fluxes(j_Fluxes), 
i_plus_inviscid_Jacobians(i_plus_inviscid_Jacobians), j_plus_inviscid_Jacobians(j_plus_inviscid_Jacobians), i_minus_inviscid_Jacobians(i_minus_inviscid_Jacobians), 
j_minus_inviscid_Jacobians(j_minus_inviscid_Jacobians), i_viscous_Jacobians(i_viscous_Jacobians), j_viscous_Jacobians(j_viscous_Jacobians),  gridtype(gridtype), dt(dt), 
CFL(CFL), inner_residual(inner_residual), outer_residual(outer_residual), cell_thermo(cell_thermo), chem_lookup_table(chem_lookup_table), progress_update(progress_update){

	V_inlet = Vector(4); U_inlet = Vector(4); 
	U = Tensor(Nx + 2, Matrix(Ny + 2, Vector(4, 0.0)));
	dU_new = Tensor(Nx, Matrix(Ny, Vector(4, 0.0)));
	dU_old = Tensor(Nx, Matrix(Ny, Vector(4, 0.0)));

	i_Fluxes = Tensor(Nx + 1, Matrix(Ny, Vector(4, 0.0)));
	j_Fluxes = Tensor(Nx, Matrix(Ny + 1, Vector(4, 0.0)));

	i_plus_inviscid_Jacobians = Tesseract(Nx + 1, Tensor(Ny, Matrix(4, Vector(4, 0.0))));
	i_minus_inviscid_Jacobians = Tesseract(Nx + 1, Tensor(Ny, Matrix(4, Vector(4, 0.0)))); 
	i_viscous_Jacobians = Tesseract(Nx + 1, Tensor(Ny, Matrix(4, Vector(4, 0.0)))); 

	j_plus_inviscid_Jacobians = Tesseract(Nx, Tensor(Ny + 1, Matrix(4, Vector(4, 0.0))));
	j_minus_inviscid_Jacobians = Tesseract(Nx, Tensor(Ny + 1, Matrix(4, Vector(4, 0.0))));
	j_viscous_Jacobians = Tesseract(Nx, Tensor(Ny + 1, Matrix(4, Vector(4, 0.0))));  

	outer_residual = 1.0;
	inner_residual = 1.0; 

	chem_lookup_table = vector<vector<ThermoEntry>>(200, vector<ThermoEntry>(200));
	cell_thermo = vector<vector<ThermoEntry>>(Nx + 2, vector<ThermoEntry>(Ny + 2));  

	V_inlet = { INLET.rho, INLET.u, INLET.v, INLET.p };

	gridtype;   
	if (dynamic_cast<RampGrid*>(&grid)) gridtype = "Ramp";  
	else if (dynamic_cast<CylinderGrid*>(&grid)) gridtype = "Cylinder";  
	else if (dynamic_cast<FlatPlateGrid*>(&grid)) gridtype = "Flat Plate"; 
	else if (dynamic_cast<DoubleConeGrid*>(&grid)) gridtype = "Double Cone"; 
	else if (dynamic_cast<MirroredGrid*>(&grid)) gridtype = "Mirrored Double Ramp";
	else gridtype = "Unknown";  

	U_inlet = primtoCons(V_inlet, perfGamma);    

	for (int i = 0; i < Nx + 2; ++i) {
		for (int j = 0; j < Ny + 2; ++j) {
			U[i][j] = U_inlet;  
		}
	}
	
	string chemistry_table = "C:/Users/Connor/source/repos/Directed Study Spring 2025/Chemical Equilibrium/Chemical_Equilibrium_Lookup_Table.csv";  // path to your CSV file
	chem_lookup_table = load_csv(chemistry_table);      

};   

Matrix Solver::inviscid_boundary_2D_E(BoundaryCondition type, const Vector& U, const Point& normals) {     

	Matrix E = zeros(4, 4);  	

	switch (type) {

	case BoundaryCondition::Inlet: 
		return E;  

	case BoundaryCondition::Outlet: 
		
		return identity(4);  

	case BoundaryCondition::IsothermalWall:  
	 
		E = { {1, 0, 0, 0}, 
			  {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
			  { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
			  { 0, 0, 0, 1 }
		};

		return E;

	case BoundaryCondition::AdiabaticWall:

		E = { {1, 0, 0, 0}, 
			  {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
			  { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
			  { 0, 0, 0, 1 }
		};

		return E;

	case BoundaryCondition::Symmetry: 		

		E = { {1, 0, 0, 0},
			  {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
			  { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
			  { 0, 0, 0, 1 }
		};

		return E; 
	default:
		throw invalid_argument("Unknown boundary condition type.");
	}

}
Vector Solver::inviscid_boundary_2D_U(BoundaryCondition type, const Vector& U, const Point& normals) {

	Vector ghost(4);
	double u = 0.0, v = 0.0;

	switch (type) {

	case BoundaryCondition::Inlet:
		return U_inlet; 

	case BoundaryCondition::Outlet:
		
		return U;

	case BoundaryCondition::IsothermalWall: 

		u = U[1] / U[0];
		v = U[2] / U[0];

		ghost[0] = U[0];
		ghost[1] = U[0] * (u - 2 * (u * normals.x + v * normals.y) * normals.x);
		ghost[2] = U[0] * (v - 2 * (u * normals.x + v * normals.y) * normals.y);
		ghost[3] = U[3];

		return ghost;

	case BoundaryCondition::AdiabaticWall:

		u = U[1] / U[0];
		v = U[2] / U[0];

		ghost[0] = U[0];
		ghost[1] = U[0] * (u - 2 * (u * normals.x + v * normals.y) * normals.x);
		ghost[2] = U[0] * (v - 2 * (u * normals.x + v * normals.y) * normals.y);
		ghost[3] = U[3];

		return ghost;

	case BoundaryCondition::Symmetry:

		u = U[1] / U[0]; 
		v = U[2] / U[0]; 

		ghost[0] = U[0]; 
		ghost[1] = U[0] * (u - 2 * (u * normals.x + v * normals.y) * normals.x);
		ghost[2] = U[0] * (v - 2 * (u * normals.x + v * normals.y) * normals.y);
		ghost[3] = U[3];

		return ghost;

	default:
		throw invalid_argument("Unknown boundary condition type.");
	}

}

void Solver::get_chemistry() {

	double e = 0;

	for (int i = 0; i < Nx + 2; ++i) { 
		for (int j = 0; j < Ny + 2; ++j) { 
			
			e = computeInternalEnergy(U[i][j]);
	
			if (e < 500000) { 
				cell_thermo[i][j].gamma = 1.4;  
				cell_thermo[i][j].R = 287;  
				cell_thermo[i][j].cv = 717;  
				cell_thermo[i][j].dpdrho = (cell_thermo[i][j].gamma - 1) * e; 
				cell_thermo[i][j].dpde = (cell_thermo[i][j].gamma - 1) * U[i][j][0];  
				cell_thermo[i][j].p = (cell_thermo[i][j].gamma - 1) * e * U[i][j][0];  
				cell_thermo[i][j].rho = U[i][j][0];  
				cell_thermo[i][j].e = e; 
				cell_thermo[i][j].T = cell_thermo[i][j].p / (cell_thermo[i][j].rho * cell_thermo[i][j].R); 
			} 
			else {  
				cell_thermo[i][j] = bilinear_interpolate(chem_lookup_table, U[i][j][0], e); 
			}		

		}
	}
}
void Solver::get_ghost_cells() {
	for (int i = 0; i < Nx; ++i) {
		U[i + 1][0] = inviscid_boundary_2D_U(BoundaryType.bottom, U[i + 1][1], grid.jNorms(i, 0));    
		U[i + 1][Ny + 1] = inviscid_boundary_2D_U(BoundaryType.top, U[i + 1][Ny], grid.jNorms(i, Ny));   
	}

	for (int j = 0; j < Ny; ++j) {
		U[0][j + 1] = inviscid_boundary_2D_U(BoundaryType.left, U[1][j + 1], grid.iNorms(0, j)); 
		U[Nx + 1][j + 1] = inviscid_boundary_2D_U(BoundaryType.right, U[Nx][j + 1], grid.iNorms(Nx, j)); 
	}
}
void Solver::compute_dt() {  

	Vector V(4); 
	double dx, dy, c, dt_old; 
	dt = 1; 

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			V = constoPrim(U[i + 1][j + 1], cell_thermo[i + 1][j + 1].gamma);  
			dx = min(grid.jArea(i, j), grid.jArea(i, j + 1));
			dy = min(grid.iArea(i, j), grid.iArea(i + 1, j));
			c = computeSoundSpeed(U[i + 1][j + 1], cell_thermo[i + 1][j + 1]); 

			dt_old = CFL / (fabs(V[1] / dx) + fabs(V[2] / dy) + c * sqrt(1 / (dx * dx) + 1 / (dy * dy))); 

			if (dt_old < dt) dt = dt_old; 

		}
	}
}

void Solver::solve_inviscid () { 

	cout << "\033[33mRunning Inviscid Chemical Equilibrium DPLR for " << Nx << " by " << Ny << " " << gridtype << " with a CFL of " << fixed << setprecision(2) << CFL << "...\033[0m" << "\n\n";
	string filename = "Inviscid " + to_string(Nx) + "x" + to_string(Ny) + "_" + gridtype + "_Solution.csv"; 

	auto start = TIME; 
	int counter = 0;
	
	t_tot = 0.0; 

	t.push_back(0.0); 
	iteration.push_back(0.0); 
	Global_Residual.push_back(10.0); 

	get_chemistry();

	while (outer_residual >= 1e-6) { 		

		compute_dt();
		
		solve_inviscid_timestep();
		compute_outer_residual();
	
		if (counter == 0) outer_residual = 1.0;		
		counter++; 

		if (counter % progress_update == 0) { 
			auto end = TIME;
			DURATION duration = end - start;
			cout << "Iteration: " << counter << "\t Residual: " << fixed << scientific << setprecision(3) << outer_residual 
				 << fixed << setprecision(2) << "\tElapsed time: " << duration.count() << " seconds" << endl;
		}

		if (counter % 1000 == 0) { 
			write_2d_csv(filename); 
		}

		t_tot += dt;
		t.push_back(t_tot);  
		Global_Residual.push_back(outer_residual); 
		iteration.push_back(counter); 

	}

	write_2d_csv(filename);  
	write_residual_csv();  

	auto end = TIME;
	DURATION duration = end - start; 
	cout << "Program completed in " << duration.count() << " seconds." << endl;


}  
void Solver::solve_inviscid_timestep() {

	inner_residual = 1.0; 

	while (inner_residual >= 1e-8) {	

		get_ghost_cells();
		get_chemistry(); 
		compute_inviscid_jacobians(); 

		solve_left_line_inviscid();

		for (int i = 1; i < Nx - 1; ++i) {
			solve_middle_line_inviscid(i);
		}

		solve_right_line_inviscid();
			
		compute_inner_residual(); 
		swap(dU_old, dU_new); 
	}

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < 4; ++k) {
				U[i + 1][j + 1][k] += dU_old[i][j][k];    
			}
		}
	}

}

void Solver::compute_inviscid_jacobians() {

	static Vector Ui(4), Uii(4), Up(4), Um(4), V1_Plus(4), V2_Plus(4), V1_Minus(4), V2_Minus(4), n(4), m(4); 
	Inviscid_State Si, Sii; 
	double g = 5.72;
	double weight, dp, pi, pii, nx, ny, l, lc, lt;
	Matrix X(4, Vector(4)), Y(4, Vector(4)); 


	// Inner i-faces  
	for (int i = 0; i <= Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {	

			nx = grid.iNorms(i, j).x;
			ny = grid.iNorms(i, j).y;

			Ui = U[i][j + 1];
			Uii = U[i + 1][j + 1];   

			pi = cell_thermo[i][j + 1].p; 
			pii = cell_thermo[i + 1][j + 1].p; 
			dp = fabs(pii - pi) / min(pii, pi);
			weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
							
			Up = weight * Ui + (1 - weight) * Uii; 
			Um = (1 - weight) * Ui + weight * Uii; 

			Si = compute_inviscid_state(Up, cell_thermo[i][j + 1], nx, ny);   
			Sii = compute_inviscid_state(Um, cell_thermo[i + 1][j + 1], nx, ny);  

			l = 0.5 * (Si.uprime + fabs(Si.uprime));
			lc = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - l);
			lt = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));

			V1_Plus = { lc * Si.k,
				(Si.u * lc + Si.a * nx * lt) * Si.k,
				(Si.v * lc + Si.a * ny * lt) * Si.k,
				(Si.h0 * lc + Si.a * Si.uprime * lt) * Si.k };

			V2_Plus = { lt / Si.a,
				Si.u * lt / Si.a + nx * lc,
				Si.v * lt / Si.a + ny * lc,
				Si.h0 * lt / Si.a + Si.uprime * lc };

			m = { Si.pp, -Si.u * Si.pe, -Si.v * Si.pe, Si.pe };   
			n = { -Si.uprime, nx, ny, 0 }; 

			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity(4);
			

			m = { Sii.pp, -Sii.u * Sii.pe, -Sii.v * Sii.pe, Sii.pe };  
			n = { -Sii.uprime, nx, ny, 0 }; 

			l = 0.5 * (Sii.uprime - fabs(Sii.uprime));
			lc = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - l);
			lt = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

			V1_Minus = { lc * Sii.k,
				(Sii.u * lc + Sii.a * nx * lt) * Sii.k,
				(Sii.v * lc + Sii.a * ny * lt) * Sii.k,
				(Sii.h0 * lc + Sii.a * Sii.uprime * lt) * Sii.k };

			V2_Minus = { lt / Sii.a,
				Sii.u * lt / Sii.a + nx * lc,
				Sii.v * lt / Sii.a + ny * lc,
				Sii.h0 * lt / Sii.a + Sii.uprime * lc };
			
			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity(4); 
			 
			i_plus_inviscid_Jacobians[i][j] = X;
			i_minus_inviscid_Jacobians[i][j] = Y;
			i_Fluxes[i][j] = X * Ui + Y * Uii;   
		}
	}


	// Inner j-faces 
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j <= Ny; ++j) {
			
			nx = grid.jNorms(i, j).x;
			ny = grid.jNorms(i, j).y;

			Ui = U[i + 1][j]; 
			Uii = U[i + 1][j + 1];   
		
			pi = cell_thermo[i + 1][j].p;  
			pii = cell_thermo[i + 1][j + 1].p; 
			dp = fabs(pii - pi) / min(pii, pi);
			weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1)); 		

			Up = weight * Ui + (1 - weight) * Uii; 
			Um = (1 - weight) * Ui + weight * Uii; 

			Si = compute_inviscid_state(Up, cell_thermo[i + 1][j], nx, ny);  
			Sii = compute_inviscid_state(Um, cell_thermo[i + 1][j + 1], nx, ny);  

			l = 0.5 * (Si.uprime + fabs(Si.uprime));
			lc = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - l); 
			lt = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a))); 

			V1_Plus = { lc * Si.k,
				(Si.u * lc + Si.a * nx * lt) * Si.k,
				(Si.v * lc + Si.a * ny * lt) * Si.k,
				(Si.h0 * lc + Si.a * Si.uprime * lt) * Si.k };

			V2_Plus = { lt / Si.a,
				Si.u * lt / Si.a + nx * lc,
				Si.v * lt / Si.a + ny * lc,
				Si.h0 * lt / Si.a + Si.uprime * lc };

			m = { Si.pp, -Si.u * Si.pe, -Si.v * Si.pe, Si.pe }; 
			n = { -Si.uprime, nx, ny, 0 };

			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + l * identity(4);

			m = { Sii.pp, -Sii.u * Sii.pe, -Sii.v * Sii.pe, Sii.pe };  
			n = { -Sii.uprime, nx, ny, 0 };

			l = 0.5 * (Sii.uprime - fabs(Sii.uprime));
			lc = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - l);
			lt = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));

			V1_Minus = { lc * Sii.k,
				(Sii.u * lc + Sii.a * nx * lt) * Sii.k,
				(Sii.v * lc + Sii.a * ny * lt) * Sii.k,
				(Sii.h0 * lc + Sii.a * Sii.uprime * lt) * Sii.k };

			V2_Minus = { lt / Sii.a,
				Sii.u * lt / Sii.a + nx * lc,
				Sii.v * lt / Sii.a + ny * lc,
				Sii.h0 * lt / Sii.a + Sii.uprime * lc };

			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + l * identity(4); 

			j_plus_inviscid_Jacobians[i][j] = X; 
			j_minus_inviscid_Jacobians[i][j] = Y;
			j_Fluxes[i][j] = X * Ui + Y * Uii;   
		}
	}

} 

void Solver::solve_left_line_inviscid() {

	int i = 0; 
	Matrix Et(4, Vector(4)), Eb(4, Vector(4)), El(4, Vector(4));
	static Matrix A(4, Vector(4)); // Relates to U(j+1) 
	static Matrix B(4, Vector(4)); // Relates to U(j) 
	static Matrix C(4, Vector(4)); // Relates to U(j-1) 
	static Vector F(4); // Right hand side  

	// Intermediate matrices
	static Matrix alpha(4, Vector(4));
	static Matrix v(Ny, Vector(4));
	static Tensor g(Ny, Matrix(4, Vector(4))); 

	// Grab boundary values
	Eb = inviscid_boundary_2D_E(BoundaryType.bottom, U[i + 1][1], grid.jNorms(i, 0));      
	Et = inviscid_boundary_2D_E(BoundaryType.top, U[i + 1][Ny], grid.jNorms(i, Ny));     
	El = inviscid_boundary_2D_E(BoundaryType.left, U[i + 1][Ny], grid.iNorms(i, Ny - 1)); 

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity(4)
		- (i_minus_inviscid_Jacobians[i][Ny - 1] + El * i_plus_inviscid_Jacobians[i][Ny - 1]) * grid.iArea(i, Ny - 1)	// Left 
		+ (i_plus_inviscid_Jacobians[i + 1][Ny - 1]) * grid.iArea(i + 1, Ny - 1)										// Right
		- (j_minus_inviscid_Jacobians[i][Ny - 1]) * grid.jArea(i, Ny - 1)												// Bottom
		+ (j_plus_inviscid_Jacobians[i][Ny] + Et * j_minus_inviscid_Jacobians[i][Ny]) * grid.jArea(i, Ny);				// Top 

	C = j_plus_inviscid_Jacobians[i][Ny - 1] * (-grid.jArea(i, Ny - 1));

	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1))
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny))
		+ i_minus_inviscid_Jacobians[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1)) * dU_old[i + 1][Ny - 1];

	alpha = A;
	v[Ny - 1] = F / alpha;
	g[Ny - 1] = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		El = inviscid_boundary_2D_E(BoundaryType.left, U[i + 1][j + 1], grid.iNorms(i, j));  

		B = j_minus_inviscid_Jacobians[i][j + 1] * (grid.jArea(i, j + 1));

		A = grid.Volume(i, j) / dt * identity(4)
			- (i_minus_inviscid_Jacobians[i][j] + El * i_plus_inviscid_Jacobians[i][j]) * grid.iArea(i, j)	// Left
			+ i_plus_inviscid_Jacobians[i + 1][j] * grid.iArea(i + 1, j)									// Right
			- j_minus_inviscid_Jacobians[i][j] * grid.jArea(i, j)											// Bottom
			+ j_plus_inviscid_Jacobians[i][j + 1] * grid.jArea(i, j + 1);									// Top

		C = j_plus_inviscid_Jacobians[i][j] * (-grid.jArea(i, j));

		F = i_Fluxes[i][j] * (grid.iArea(i, j))
			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j))
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1))
			+ i_minus_inviscid_Jacobians[i + 1][j] * (-grid.iArea(i + 1, j)) * dU_old[i + 1][j];



		alpha = A - B * g[j + 1];
		g[j] = C / alpha;
		v[j] = (F - B * v[j + 1]) / alpha;
	}

	//  Bottom boundary

	El = inviscid_boundary_2D_E(BoundaryType.left, U[i + 1][1], grid.iNorms(i, 0));   

	B = j_minus_inviscid_Jacobians[i][1] * (grid.jArea(i, 1));

	A = grid.Volume(i, 0) / dt * identity(4)
		- (i_minus_inviscid_Jacobians[i][0] + El * i_plus_inviscid_Jacobians[i][0]) * grid.iArea(i, 0)		// Left
		+ i_plus_inviscid_Jacobians[i + 1][0] * grid.iArea(i + 1, 0)										// Right
		- (Eb * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0]) * grid.jArea(i, 0)		// Bottom
		+ j_plus_inviscid_Jacobians[i][1] * grid.jArea(i, 1);												// Top

	F = i_Fluxes[i][0] * (grid.iArea(i, 0))
		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1))
		+ i_minus_inviscid_Jacobians[i + 1][0] * (-grid.iArea(i + 1, 0)) * dU_old[i + 1][0];


	alpha = A - B * g[1];
	v[0] = (F - B * v[1]) / alpha;

	// Calculate dU
	dU_new[i][0] = v[0];
	for (int j = 1; j < Ny; ++j) {
		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];
	}

}

void Solver::solve_middle_line_inviscid(const int i) {

	//auto start = TIME;
	Matrix Et(4, Vector(4)), Eb(4, Vector(4));
	static Matrix A(4, Vector(4)); // Relates to U(j+1) 
	static Matrix B(4, Vector(4)); // Relates to U(j) 
	static Matrix C(4, Vector(4)); // Relates to U(j-1) 
	static Vector F(4); // Right hand side  

	// Intermediate matrices
	static Matrix alpha(4, Vector(4));
	static Matrix v(Ny, Vector(4));
	static Tensor g(Ny, Matrix(4, Vector(4)));

	// Grab boundary values
	Eb = inviscid_boundary_2D_E(BoundaryType.bottom, U[i + 1][1], grid.jNorms(i, 0));  
	Et = inviscid_boundary_2D_E(BoundaryType.top, U[i + 1][Ny], grid.jNorms(i, Ny));    

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity(4) 
		- i_minus_inviscid_Jacobians[i][Ny - 1] * grid.iArea(i, Ny - 1) 
		+ i_plus_inviscid_Jacobians[i + 1][Ny - 1] * grid.iArea(i + 1, Ny - 1)
		- j_minus_inviscid_Jacobians[i][Ny - 1] * grid.jArea(i, Ny - 1) 
		+ (j_plus_inviscid_Jacobians[i][Ny] + Et * j_minus_inviscid_Jacobians[i][Ny]) * grid.jArea(i, Ny); 

	C = j_plus_inviscid_Jacobians[i][Ny - 1] * (-grid.jArea(i, Ny - 1)); 

	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1)) 
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny))
		+ i_plus_inviscid_Jacobians[i][Ny - 1] * (grid.iArea(i, Ny - 1)) * dU_old[i - 1][Ny - 1]
		+ i_minus_inviscid_Jacobians[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1)) * dU_old[i + 1][Ny - 1]; 
	
	alpha = A;
	v[Ny - 1] = F / alpha;  
	g[Ny - 1] = C / alpha; 

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		B = j_minus_inviscid_Jacobians[i][j + 1] * (grid.jArea(i, j + 1)); 

		A = grid.Volume(i, j) / dt * identity(4) 
			- i_minus_inviscid_Jacobians[i][j] * grid.iArea(i, j) 
			+ i_plus_inviscid_Jacobians[i + 1][j] * grid.iArea(i + 1, j)  
			- j_minus_inviscid_Jacobians[i][j] * grid.jArea(i, j)  
			+ j_plus_inviscid_Jacobians[i][j + 1] * grid.jArea(i, j + 1);   

		C = j_plus_inviscid_Jacobians[i][j] * (-grid.jArea(i, j));  

		F = i_Fluxes[i][j] * (grid.iArea(i, j))  
			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j)) 
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1)) 
			+ i_plus_inviscid_Jacobians[i][j] * (grid.iArea(i, j)) * dU_old[i - 1][j] 
			+ i_minus_inviscid_Jacobians[i + 1][j] * (-grid.iArea(i + 1, j)) * dU_old[i + 1][j]; 



		alpha = A - B * g[j + 1]; 
		g[j] = C / alpha; 
		v[j] = (F - B * v[j + 1]) / alpha; 
	}

	//  Bottom boundary
	B = j_minus_inviscid_Jacobians[i][1] * (grid.jArea(i, 1));

	A = grid.Volume(i, 0) / dt * identity(4)
		- i_minus_inviscid_Jacobians[i][0] * grid.iArea(i, 0)
		+ i_plus_inviscid_Jacobians[i + 1][0] * grid.iArea(i + 1, 0)
		- (Eb * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0]) * grid.jArea(i, 0) 
		+ j_plus_inviscid_Jacobians[i][1] * grid.jArea(i, 1); 

	F = i_Fluxes[i][0] * (grid.iArea(i, 0))
		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1))
		+ i_plus_inviscid_Jacobians[i][0] * (grid.iArea(i, 0)) * dU_old[i - 1][0]
		+ i_minus_inviscid_Jacobians[i + 1][0] * (-grid.iArea(i + 1, 0)) * dU_old[i + 1][0];


	alpha = A - B * g[1];
	v[0] = (F - B * v[1]) / alpha; 

	// Calculate dU
	dU_new[i][0] = v[0]; 

	for (int j = 1; j < Ny; ++j) {
		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];  
	}

	//auto end = TIME;
	//DURATION duration = end - start; 
	//cout << duration.count() << endl; 
}

void Solver::solve_right_line_inviscid() {

	int i = Nx - 1; 
	Matrix Et(4, Vector(4)), Eb(4, Vector(4)), Er(4, Vector(4));

	static Matrix A(4, Vector(4)); // Relates to U(j+1) 
	static Matrix B(4, Vector(4)); // Relates to U(j) 
	static Matrix C(4, Vector(4)); // Relates to U(j-1) 
	static Vector F(4); // Right hand side  

	// Intermediate matrices
	static Matrix alpha(4, Vector(4));
	static Matrix v(Ny, Vector(4));
	static Tensor g(Ny, Matrix(4, Vector(4)));

	// Grab boundary values
	Eb = inviscid_boundary_2D_E(BoundaryType.bottom, U[i + 1][1], grid.jNorms(i, 0));
	Et = inviscid_boundary_2D_E(BoundaryType.top, U[i + 1][Ny], grid.jNorms(i, Ny));
	Er = inviscid_boundary_2D_E(BoundaryType.right, U[i + 1][Ny], grid.iNorms(i + 1, Ny - 1)); 

	// Top Boundary
	A = grid.Volume(i, Ny - 1) / dt * identity(4)
		- i_minus_inviscid_Jacobians[i][Ny - 1] * grid.iArea(i, Ny - 1)																// Left
		+ (i_plus_inviscid_Jacobians[i + 1][Ny - 1] + Er * i_minus_inviscid_Jacobians[i + 1][Ny - 1]) * grid.iArea(i + 1, Ny - 1)	// Right
		- j_minus_inviscid_Jacobians[i][Ny - 1] * grid.jArea(i, Ny - 1)																// Bottom
		+ (j_plus_inviscid_Jacobians[i][Ny] + Et * j_minus_inviscid_Jacobians[i][Ny]) * grid.jArea(i, Ny);							// Top 

	C = j_plus_inviscid_Jacobians[i][Ny - 1] * (-grid.jArea(i, Ny - 1));

	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1))
		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny))
		+ i_plus_inviscid_Jacobians[i][Ny - 1] * (grid.iArea(i, Ny - 1)) * dU_old[i - 1][Ny - 1];

	alpha = A;
	v[Ny - 1] = F / alpha;
	g[Ny - 1] = C / alpha;

	// Middle Boundry
	for (int j = Ny - 2; j > 0; --j) {

		Er = inviscid_boundary_2D_E(BoundaryType.right, U[i + 1][j + 1], grid.iNorms(i + 1, j));   

		B = j_minus_inviscid_Jacobians[i][j + 1] * (grid.jArea(i, j + 1));

		A = grid.Volume(i, j) / dt * identity(4)
			- i_minus_inviscid_Jacobians[i][j] * grid.iArea(i, j)														// Left
			+ (i_plus_inviscid_Jacobians[i + 1][j] + Er * i_minus_inviscid_Jacobians[i + 1][j])* grid.iArea(i + 1, j)	// Right
			- j_minus_inviscid_Jacobians[i][j] * grid.jArea(i, j)														// Bottom
			+ j_plus_inviscid_Jacobians[i][j + 1] * grid.jArea(i, j + 1);												// Top

		C = j_plus_inviscid_Jacobians[i][j] * (-grid.jArea(i, j));

		F = i_Fluxes[i][j] * (grid.iArea(i, j))
			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
			+ j_Fluxes[i][j] * (grid.jArea(i, j))
			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1))
			+ i_plus_inviscid_Jacobians[i][j] * (grid.iArea(i, j)) * dU_old[i - 1][j];



		alpha = A - B * g[j + 1];
		g[j] = C / alpha;
		v[j] = (F - B * v[j + 1]) / alpha;
	}
	//  Bottom boundary
	Er = inviscid_boundary_2D_E(BoundaryType.right, U[i + 1][1], grid.iNorms(i + 1, 0)); 

	B = j_minus_inviscid_Jacobians[i][1] * (grid.jArea(i, 1));

	A = grid.Volume(i, 0) / dt * identity(4)
		- i_minus_inviscid_Jacobians[i][0] * grid.iArea(i, 0)															// Left
		+ (i_plus_inviscid_Jacobians[i + 1][0] + Er * i_minus_inviscid_Jacobians[i + 1][0]) * grid.iArea(i + 1, 0)		// Right
		- (Eb * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0]) * grid.jArea(i, 0)					// Bottom
		+ j_plus_inviscid_Jacobians[i][1] * grid.jArea(i, 1);															// Top 

	F = i_Fluxes[i][0] * (grid.iArea(i, 0))
		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0))
		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
		+ j_Fluxes[i][1] * (-grid.jArea(i, 1))
		+ i_plus_inviscid_Jacobians[i][0] * (grid.iArea(i, 0)) * dU_old[i - 1][0];


	alpha = A - B * g[1];
	v[0] = (F - B * v[1]) / alpha;

	// Calculate dU
	dU_new[i][0] = v[0];

	for (int j = 1; j < Ny; ++j) {
		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];
	}

	

}

void Solver::compute_inner_residual() { 
	
	inner_residual = 0.0;
	double res = 0.0; 

	static Vector A(4); // Relates to U(j+1) 
	static Vector B(4); // Relates to U(j) 
	static Vector C(4); // Relates to U(j-1) 
	static double F; // Right hand side  

	Vector one = { 1, 0, 0, 0 };

	// Middle Boundary
	for (int i = 1; i < Nx - 1; ++i) {
		for (int j = 1; j < Ny - 1; ++j) {

			B = j_minus_inviscid_Jacobians[i][j + 1][0] * (grid.jArea(i, j + 1)); 

			A = grid.Volume(i, j) / dt * one 
				- i_minus_inviscid_Jacobians[i][j][0] * grid.iArea(i, j) 
				+ i_plus_inviscid_Jacobians[i + 1][j][0] * grid.iArea(i + 1, j)
				- j_minus_inviscid_Jacobians[i][j][0] * grid.jArea(i, j) 
				+ j_plus_inviscid_Jacobians[i][j + 1][0] * grid.jArea(i, j + 1);

			C = j_plus_inviscid_Jacobians[i][j][0] * (-grid.jArea(i, j)); 

			F = i_Fluxes[i][j][0] * (grid.iArea(i, j)) 
				+ i_Fluxes[i + 1][j][0] * (-grid.iArea(i + 1, j)) 
				+ j_Fluxes[i][j][0] * (grid.jArea(i, j)) 
				+ j_Fluxes[i][j + 1][0] * (-grid.jArea(i, j + 1)) 
				+ i_plus_inviscid_Jacobians[i][j][0] * (grid.iArea(i, j)) * dU_old[i - 1][j]
				+ i_minus_inviscid_Jacobians[i + 1][j][0] * (-grid.iArea(i + 1, j)) * dU_old[i + 1][j];   

			res = B * dU_new[i][j + 1] + A * dU_new[i][j] + C * dU_new[i][j - 1] - F;
			inner_residual += res * res;  

		}
	}
	inner_residual = sqrt(inner_residual); 
}

void Solver::compute_outer_residual() {  

	outer_residual = 0.0;  
	double intres;

	for (int i = 1; i < Nx - 1; ++i) {
		for (int j = 1; j < Ny - 1; ++j) {
			intres = (i_Fluxes[i][j][0] * (-grid.iArea(i, j))			// Left Side		   
				+ i_Fluxes[i + 1][j][0] * (grid.iArea(i + 1, j))		// Right Side 
				+ j_Fluxes[i][j][0] * (-grid.jArea(i, j))			// Bottom Side  
				+ j_Fluxes[i][j + 1][0] * (grid.jArea(i, j + 1)))/ grid.Volume(i, j); 

			outer_residual = outer_residual + intres * intres; 
		
		}
	}

	outer_residual = sqrt(outer_residual);
}

Vector Solver::minmod(Vector& Ui, Vector& Uii) {
	
	Vector result(4);  

	for (int i = 0; i < 4; ++i) {
		if (Ui[i] * Uii[i] <= 0) {
			result[i] = 0;
		}
		else result[i] = min(Ui[i], Uii[i]);
	}

	return result; 

}

void Solver::write_residual_csv() {

	string filename = "RESIDUAL_PLOT.csv";

	ofstream file(filename);
	file << "Residual, time, iteration" << endl;
	for (int i = 0; i < Global_Residual.size(); ++i) {
		file << Global_Residual[i] << ", " << t[i] << ", " << iteration[i] << endl;
	}
	file.close();
	cout << "\033[36mResidual File saved successfully as \"" << filename << "\"\033[0m" << endl;

}

void Solver::write_2d_csv(const string& filename) {

	double a, density, pressure, u_velocity, v_velocity, Mach, Temperature;

	Vector Primitives(4);

	ofstream file(filename);
	file << "density, u velocity, v velocity, pressure, Mach, Temperature, x_points, y_points, z_points" << endl;

	get_chemistry();  
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {

			Primitives = constoPrim(U[i + 1][j + 1], cell_thermo[i + 1][j + 1].gamma); 
			density = Primitives[0];
			pressure = Primitives[3];
			a = computeSoundSpeed(U[i + 1][j + 1], cell_thermo[i + 1][j + 1]); 
			u_velocity = Primitives[1];
			v_velocity = Primitives[2];
			Mach = sqrt(u_velocity * u_velocity + v_velocity * v_velocity) / a;
			Temperature = pressure / (density * cell_thermo[i + 1][j + 1].R); 

			file << density << ", " << u_velocity << ", " << v_velocity << ", " << pressure << ", " << Mach << ", " << Temperature<< ", " << grid.Center(i, j).x << ", " << grid.Center(i, j).y << ", 0.0" << endl;

		}
	}

	file.close();
	cout << "\033[36m2D File saved successfully as \"" << filename << "\"\033[0m" << endl;  
}



//
//Vector Solver::constoViscPrim(const Vector& U, ThermoEntry& Thermo) {   
//Vector result(4);
//result[0] = U[0];
//result[1] = U[1] / U[0];
//result[2] = U[2] / U[0];
//result[3] = computeTemperature(U, Thermo);
//return result;
//}
// 
// Matrix Solver::viscous_boundary_2D_E(BoundaryCondition type, const Vector& U, const Point& normals) {
//
//	Matrix E = zeros(4, 4);
//	double Ti;
//
//	switch (type) {
//
//	case BoundaryCondition::Inlet:
//		return E;
//
//	case BoundaryCondition::Outlet:
//
//		return identity(4);
//
//	case BoundaryCondition::IsothermalWall:
//
//		Ti = computeTemperature(U);
//
//		E = { {2 * Ti / Tw - 1 , 0, 0, 2 * U[0] / Tw},
//			  { 0, -1, 0, 0 },
//			  { 0, 0, -1, 0},
//			  { 0, 0, 0, -1 }
//		};
//
//		return E;
//
//	case BoundaryCondition::AdiabaticWall:
//
//		Ti = computeTemperature(U);
//
//		E = { {2 * Ti / Tw - 1 , 0, 0, 2 * U[0] / Tw},
//			  { 0, -1, 0, 0 },
//			  { 0, 0, -1, 0},
//			  { 0, 0, 0, 1 }
//		};
//
//		return E;
//
//	case BoundaryCondition::Symmetry:
//
//		E = { {1, 0, 0, 0},
//			  {0, (1 - 2 * normals.x * normals.x), -2 * normals.x * normals.y, 0 },
//			  { 0, -2 * normals.x * normals.y, (1 - 2 * normals.y * normals.y), 0 },
//			  { 0, 0, 0, 1 }
//		};
//
//		return E;
//	default:
//		throw invalid_argument("Unknown boundary condition type.");
//	}
//
//}
//Vector Solver::viscous_boundary_2D_U(BoundaryCondition type, const Vector& U, const Point& normals) {
//
//	Vector ghost(4);
//	double T_in;
//	double u = 0.0, v = 0.0;
//
//	switch (type) {
//
//	case BoundaryCondition::Inlet:
//		return U_inlet;
//
//	case BoundaryCondition::Outlet:
//
//		return U;
//
//	case BoundaryCondition::IsothermalWall:
//
//		T_in = computeTemperature(U);
//
//		ghost[0] = 2 * (U[0] * T_in / Tw) - U[0];
//		ghost[1] = -U[1];
//		ghost[2] = -U[2];
//		ghost[3] = U[3];
//
//		return ghost;
//
//	case BoundaryCondition::AdiabaticWall:
//
//		T_in = Tw;
//
//		ghost[0] = 2 * (U[0] * T_in / Tw) - U[0];
//		ghost[1] = -U[1];
//		ghost[2] = -U[2];
//		ghost[3] = U[3];
//
//		return ghost;
//
//	case BoundaryCondition::Symmetry:
//
//		u = U[1] / U[0];
//		v = U[2] / U[0];
//
//		ghost[0] = U[0];
//		ghost[1] = U[0] * (u - 2 * (u * normals.x + v * normals.y) * normals.x);
//		ghost[2] = U[0] * (v - 2 * (u * normals.x + v * normals.y) * normals.y);
//		ghost[3] = U[3];
//
//		return ghost;
//
//	default:
//		throw invalid_argument("Unknown boundary condition type.");
//	}
//
//}
//void Solver::solve_viscous() {
//
//	cout << "\033[33mRunning Viscous DPLR for " << Nx << " by " << Ny << " " << gridtype << " with a CFL of " << fixed << setprecision(2) << CFL << "...\033[0m" << "\n\n";
//	string filename = "Viscous " + to_string(Nx) + "x" + to_string(Ny) + "_" + gridtype + "_Solution.csv";
//
//	auto start = TIME;
//	int counter = 0;
//
//	t_tot = 0.0;
//
//	t.push_back(0.0);
//	iteration.push_back(0.0);
//	Global_Residual.push_back(0.0);
//
//	while (outer_residual >= 1e-6) {
//
//		compute_dt();
//
//		solve_viscous_timestep();
//
//		compute_outer_residual();
//
//		if (counter == 0) outer_residual = 1.0;
//		counter++;
//
//
//		if (counter % progress_update == 0) {
//			auto end = TIME;
//			DURATION duration = end - start;
//			cout << "Iteration: " << counter << "\t Residual: " << fixed << scientific << setprecision(3) << outer_residual
//				<< fixed << setprecision(2) << "\tElapsed time: " << duration.count() << " seconds" << endl;
//		}
//
//		if (counter % 1000 == 0) {
//			write_2d_csv(filename);
//		}
//
//		t_tot += dt;
//		t.push_back(t_tot);
//		Global_Residual.push_back(outer_residual);
//		iteration.push_back(counter);
//
//	}
//
//	write_2d_csv(filename);
//	write_1d_csv("1D VISCOUS.csv");
//	write_residual_csv();
//	viscous_calculations();
//
//	auto end = TIME;
//	DURATION duration = end - start;
//	cout << "Program complete in " << duration.count() << " seconds." << endl;
//}
//void Solver::solve_viscous_timestep() {
//
//	inner_residual = 1.0;
//
//	while (inner_residual >= 1e-8) {
//
//		compute_viscous_jacobians();
//
//		solve_left_line_viscous();
//
//		for (int i = 1; i < Nx - 1; ++i) {
//			solve_middle_line_viscous(i);
//		}
//
//		solve_right_line_viscous();
//
//		compute_inner_residual();
//		swap(dU_old, dU_new);
//	}
//
//	for (int i = 0; i < Nx; ++i) {
//		for (int j = 0; j < Ny; ++j) {
//			for (int k = 0; k < 4; ++k) {
//				U[i][j][k] += dU_old[i][j][k];
//			}
//		}
//	}
//}
//void Solver::compute_viscous_jacobians() {
//
//	static Vector Ui(4), Uii(4), Up(4), Um(4), V1_Plus(4), V2_Plus(4), V1_Minus(4), V2_Minus(4), n(4), m(4), dv(4), W(4);
//	Viscous_State Si, Sii;
//	double g = 5.72;
//	double weight, dp, pi, pii, pe, nx, ny, lp, lcp, ltp, lm, lcm, ltm, rho, u, v, T, mu, lambda, k, dn, dx, dy;
//	pe = (gamma - 1);
//	Matrix M(4, Vector(4)), N(4, Vector(4)), X(4, Vector(4)), Y(4, Vector(4)), Elv(4, Vector(4)), Erv(4, Vector(4)), Ebv(4, Vector(4)), Etv(4, Vector(4));
//
//	// Calculate Jacobians and Explicit fluxes for i-faces on left boundary 
//	for (int j = 0; j < Ny; ++j) {
//
//		Uii = U[0][j];
//
//		Elv = viscous_boundary_2D_E(BoundaryType.left, Uii, grid.iNorms(0, j));
//		Ui = viscous_boundary_2D_U(BoundaryType.left, Uii, grid.iNorms(0, j));
//
//		nx = grid.iNorms(0, j).x;
//		ny = grid.iNorms(0, j).y;
//
//		pi = computePressure(Ui);
//		pii = computePressure(Uii);
//		dp = fabs(pii - pi) / min(pii, pi);
//		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
//
//		Up = weight * Ui + (1 - weight) * Uii;
//		Um = (1 - weight) * Ui + weight * Uii;
//
//		Si = compute_viscous_state(Up, nx, ny);
//		Sii = compute_viscous_state(Um, nx, ny);
//
//		W = 0.5 * (Ui + Uii);
//		rho = W[0];
//		u = W[1] / W[0];
//		v = W[2] / W[0];
//		T = computeTemperature(W);
//		mu = 1.458 * 1e-6 * T * sqrt(T) / (T + 110.3);
//		lambda = -2 / 3 * mu;
//		k = cp * mu / Pr;
//
//		M[0] = { 0, 0, 0, 0 };
//		M[1][0] = 0;
//		M[1][1] = (lambda + 2 * mu) * nx * nx + mu * ny * ny;
//		M[1][2] = (lambda + mu) * nx * ny;
//		M[1][3] = 0;
//		M[2][0] = 0;
//		M[2][1] = (mu + lambda) * nx * ny;
//		M[2][2] = mu * nx * nx + (lambda + 2 * mu) * ny * ny;
//		M[2][3] = 0;
//		M[3][0] = 0;
//		M[3][1] = u * (lambda + 2 * mu) * nx * nx + (v * mu + v * lambda) * nx * ny + u * mu * ny * ny;
//		M[3][2] = v * mu * nx * nx + (u * lambda + u * mu) * nx * ny + v * (lambda + 2 * mu) * ny * ny;
//		M[3][3] = k * (nx * nx + ny * ny);
//
//
//		N = {
//			{ {1, 0, 0, 0},
//			{-u / rho, 1 / rho, 0, 0},
//			{-v / rho, 0, 1 / rho, 0},
//			{ (u * u + v * v) / (2 * cv * rho) - T / rho, -u / (cv * rho), -v / (cv * rho), 1 / (cv * rho)}
//			}
//		};
//
//		dx = grid.Center(1, j).x - grid.Center(0, j).x;
//		dy = grid.Center(1, j).y - grid.Center(0, j).y;
//		dn = sqrt(dx * dx + dy * dy);
//
//		dv = constoViscPrim(Uii) - constoViscPrim(Ui);
//
//		lp = 0.5 * (Si.uprime + fabs(Si.uprime));
//		lcp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - lp);
//		ltp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));
//
//		V1_Plus = { lcp * Si.k,
//			(Si.u * lcp + Si.a * nx * ltp) * Si.k,
//			(Si.v * lcp + Si.a * ny * ltp) * Si.k,
//			(Si.h0 * lcp + Si.a * Si.uprime * ltp) * Si.k };
//
//		V2_Plus = { ltp / Si.a,
//			Si.u * ltp / Si.a + nx * lcp,
//			Si.v * ltp / Si.a + ny * lcp,
//			Si.h0 * ltp / Si.a + Si.uprime * lcp };
//
//		m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
//		n = { -Si.uprime, nx, ny, 0 };
//
//		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(4);
//
//		m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
//		n = { -Sii.uprime, nx, ny, 0 };
//
//		lm = 0.5 * (Sii.uprime - fabs(Sii.uprime));
//		lcm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - lm);
//		ltm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));
//
//		V1_Minus = { lcm * Sii.k,
//				(Sii.u * lcm + Sii.a * nx * ltm) * Sii.k,
//				(Sii.v * lcm + Sii.a * ny * ltm) * Sii.k,
//				(Sii.h0 * lcm + Sii.a * Sii.uprime * ltm) * Sii.k };
//
//		V2_Minus = { ltm / Sii.a,
//			Sii.u * ltm / Sii.a + nx * lcm,
//			Sii.v * ltm / Sii.a + ny * lcm,
//			Sii.h0 * ltm / Sii.a + Sii.uprime * lcm };
//
//
//		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(4);
//
//		i_viscous_Jacobians[0][j] = M / (-dn) * (identity(4) - Elv) * N;
//
//		//displayMatrix(Z);  
//		i_plus_inviscid_Jacobians[0][j] = X;
//		i_minus_inviscid_Jacobians[0][j] = Y;
//		i_Fluxes[0][j] = X * Ui + Y * Uii - M * dv / dn;
//	}
//
//	// Calculate Jacobians and Explicit fluxes for i-faces on right boundary
//	for (int j = 0; j < Ny; ++j) {
//
//		Ui = U[Nx - 1][j];
//		Uii = viscous_boundary_2D_U(BoundaryType.right, Ui, grid.iNorms(Nx, j));
//		Erv = viscous_boundary_2D_E(BoundaryType.right, Ui, grid.iNorms(Nx, j));
//
//		nx = grid.iNorms(Nx, j).x;
//		ny = grid.iNorms(Nx, j).y;
//
//		pi = computePressure(Ui);
//		pii = computePressure(Uii);
//		dp = fabs(pii - pi) / min(pii, pi);
//		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
//
//		Up = weight * Ui + (1 - weight) * Uii;
//		Um = (1 - weight) * Ui + weight * Uii;
//
//		Si = compute_viscous_state(Up, nx, ny);
//		Sii = compute_viscous_state(Um, nx, ny);
//
//		W = 0.5 * (Ui + Uii);
//		rho = W[0];
//		u = W[1] / W[0];
//		v = W[2] / W[0];
//		T = computeTemperature(W);
//		mu = 1.458 * 1e-6 * T * sqrt(T) / (T + 110.3);
//		lambda = -2 / 3 * mu;
//		k = cp * mu / Pr;
//
//		M[0] = { 0, 0, 0, 0 };
//		M[1][0] = 0;
//		M[1][1] = (lambda + 2 * mu) * nx * nx + mu * ny * ny;
//		M[1][2] = (lambda + mu) * nx * ny;
//		M[1][3] = 0;
//		M[2][0] = 0;
//		M[2][1] = (mu + lambda) * nx * ny;
//		M[2][2] = mu * nx * nx + (lambda + 2 * mu) * ny * ny;
//		M[2][3] = 0;
//		M[3][0] = 0;
//		M[3][1] = u * (lambda + 2 * mu) * nx * nx + (v * mu + v * lambda) * nx * ny + u * mu * ny * ny;
//		M[3][2] = v * mu * nx * nx + (u * lambda + u * mu) * nx * ny + v * (lambda + 2 * mu) * ny * ny;
//		M[3][3] = k * (nx * nx + ny * ny);
//
//		N = {
//			{ {1, 0, 0, 0},
//			{-u / rho, 1 / rho, 0, 0},
//			{-v / rho, 0, 1 / rho, 0},
//			{ (u * u + v * v) / (2 * cv * rho) - T / rho, -u / (cv * rho), -v / (cv * rho), 1 / (cv * rho)}
//			}
//		};
//
//		dx = grid.Center(Nx - 1, j).x - grid.Center(Nx - 2, j).x;
//		dy = grid.Center(Nx - 1, j).y - grid.Center(Nx - 2, j).y;
//
//		dn = sqrt(dx * dx + dy * dy);
//		dv = constoViscPrim(Uii) - constoViscPrim(Ui);
//
//
//		lp = 0.5 * (Si.uprime + fabs(Si.uprime));
//		lcp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - lp);
//		ltp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));
//
//		V1_Plus = { lcp * Si.k,
//			(Si.u * lcp + Si.a * nx * ltp) * Si.k,
//			(Si.v * lcp + Si.a * ny * ltp) * Si.k,
//			(Si.h0 * lcp + Si.a * Si.uprime * ltp) * Si.k };
//
//		V2_Plus = { ltp / Si.a,
//			Si.u * ltp / Si.a + nx * lcp,
//			Si.v * ltp / Si.a + ny * lcp,
//			Si.h0 * ltp / Si.a + Si.uprime * lcp };
//
//		m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
//		n = { -Si.uprime, nx, ny, 0 };
//
//		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(4);
//
//		m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
//		n = { -Sii.uprime, nx, ny, 0 };
//
//		lm = 0.5 * (Sii.uprime - fabs(Sii.uprime));
//		lcm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - lm);
//		ltm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));
//
//		V1_Minus = { lcm * Sii.k,
//				(Sii.u * lcm + Sii.a * nx * ltm) * Sii.k,
//				(Sii.v * lcm + Sii.a * ny * ltm) * Sii.k,
//				(Sii.h0 * lcm + Sii.a * Sii.uprime * ltm) * Sii.k };
//
//		V2_Minus = { ltm / Sii.a,
//			Sii.u * ltm / Sii.a + nx * lcm,
//			Sii.v * ltm / Sii.a + ny * lcm,
//			Sii.h0 * ltm / Sii.a + Sii.uprime * lcm };
//
//		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(4);
//
//		i_viscous_Jacobians[Nx][j] = M / (-dn) * (Erv - identity(4)) * N;
//		i_plus_inviscid_Jacobians[Nx][j] = X;
//		i_minus_inviscid_Jacobians[Nx][j] = Y;
//		i_Fluxes[Nx][j] = X * Ui + Y * Uii - M * dv / dn;
//	}
//
//	// Calculate Jacobians and Explicit fluxes for j-faces on bottom boundary
//	for (int i = 0; i < Nx; ++i) {
//
//		Uii = U[i][0];
//		Ui = viscous_boundary_2D_U(BoundaryType.bottom, Uii, grid.jNorms(i, 0));
//		Ebv = viscous_boundary_2D_E(BoundaryType.bottom, Uii, grid.jNorms(i, 0));
//
//		nx = grid.jNorms(i, 0).x;
//		ny = grid.jNorms(i, 0).y;
//
//		pi = computePressure(Ui);
//		pii = computePressure(Uii);
//		dp = fabs(pii - pi) / min(pii, pi);
//		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
//
//		Up = weight * Ui + (1 - weight) * Uii;
//		Um = (1 - weight) * Ui + weight * Uii;
//		Si = compute_viscous_state(Up, nx, ny);
//		Sii = compute_viscous_state(Um, nx, ny);
//
//		W = 0.5 * (Ui + Uii);
//		rho = W[0];
//		u = W[1] / W[0];
//		v = W[2] / W[0];
//		T = computeTemperature(W);
//
//		mu = 1.458 * 1e-6 * T * sqrt(T) / (T + 110.3);
//		lambda = -2 / 3 * mu;
//		k = cp * mu / Pr;
//
//		M[0] = { 0, 0, 0, 0 };
//		M[1][0] = 0;
//		M[1][1] = (lambda + 2 * mu) * nx * nx + mu * ny * ny;
//		M[1][2] = (lambda + mu) * nx * ny;
//		M[1][3] = 0;
//		M[2][0] = 0;
//		M[2][1] = (mu + lambda) * nx * ny;
//		M[2][2] = mu * nx * nx + (lambda + 2 * mu) * ny * ny;
//		M[2][3] = 0;
//		M[3][0] = 0;
//		M[3][1] = u * (lambda + 2 * mu) * nx * nx + (v * mu + v * lambda) * nx * ny + u * mu * ny * ny;
//		M[3][2] = v * mu * nx * nx + (u * lambda + u * mu) * nx * ny + v * (lambda + 2 * mu) * ny * ny;
//		M[3][3] = k * (nx * nx + ny * ny);
//
//		N = {
//			{ {1, 0, 0, 0},
//			{-u / rho, 1 / rho, 0, 0},
//			{-v / rho, 0, 1 / rho, 0},
//			{ (u * u + v * v) / (2 * cv * rho) - T / rho, -u / (cv * rho), -v / (cv * rho), 1 / (cv * rho)}
//			}
//		};
//
//		dx = grid.Center(i, 1).x - grid.Center(i, 0).x;
//		dy = grid.Center(i, 1).y - grid.Center(i, 0).y;
//
//		dn = sqrt(dx * dx + dy * dy);
//		dv = constoViscPrim(Uii) - constoViscPrim(Ui);
//
//		lp = 0.5 * (Si.uprime + fabs(Si.uprime));
//		lcp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - lp);
//		ltp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));
//
//		V1_Plus = { lcp * Si.k,
//				(Si.u * lcp + Si.a * nx * ltp) * Si.k,
//				(Si.v * lcp + Si.a * ny * ltp) * Si.k,
//				(Si.h0 * lcp + Si.a * Si.uprime * ltp) * Si.k };
//
//		V2_Plus = { ltp / Si.a,
//			Si.u * ltp / Si.a + nx * lcp,
//			Si.v * ltp / Si.a + ny * lcp,
//			Si.h0 * ltp / Si.a + Si.uprime * lcp };
//
//		m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
//		n = { -Si.uprime, nx, ny, 0 };
//
//		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(4);
//
//		m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
//		n = { -Sii.uprime, nx, ny, 0 };
//
//		lm = 0.5 * (Sii.uprime - fabs(Sii.uprime));
//		lcm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - lm);
//		ltm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));
//
//		V1_Minus = { lcm * Sii.k,
//				(Sii.u * lcm + Sii.a * nx * ltm) * Sii.k,
//				(Sii.v * lcm + Sii.a * ny * ltm) * Sii.k,
//				(Sii.h0 * lcm + Sii.a * Sii.uprime * ltm) * Sii.k };
//
//		V2_Minus = { ltm / Sii.a,
//			Sii.u * ltm / Sii.a + nx * lcm,
//			Sii.v * ltm / Sii.a + ny * lcm,
//			Sii.h0 * ltm / Sii.a + Sii.uprime * lcm };
//
//		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(4);
//
//		j_viscous_Jacobians[i][0] = M / (-dn) * (identity(4) - Ebv) * N;
//		j_plus_inviscid_Jacobians[i][0] = X;
//		j_minus_inviscid_Jacobians[i][0] = Y;
//		j_Fluxes[i][0] = X * Ui + Y * Uii - M * dv / dn;
//
//	}
//
//	// Calculate Jacobians and Explicit fluxes for j-faces on top boundary
//	for (int i = 0; i < Nx; ++i) {
//
//		Ui = U[i][Ny - 1];
//		Uii = viscous_boundary_2D_U(BoundaryType.top, Ui, grid.jNorms(i, Ny));
//		Etv = viscous_boundary_2D_E(BoundaryType.top, Ui, grid.jNorms(i, Ny));
//
//		nx = grid.jNorms(i, Ny).x;
//		ny = grid.jNorms(i, Ny).y;
//
//		pi = computePressure(Ui);
//		pii = computePressure(Uii);
//		dp = fabs(pii - pi) / min(pii, pi);
//		weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
//
//		Up = weight * Ui + (1 - weight) * Uii;
//		Um = (1 - weight) * Ui + weight * Uii;
//
//		Si = compute_viscous_state(Up, nx, ny);
//		Sii = compute_viscous_state(Um, nx, ny);
//
//		W = 0.5 * (Ui + Uii);
//		rho = W[0];
//		u = W[1] / W[0];
//		v = W[2] / W[0];
//		T = computeTemperature(W);
//		mu = 1.458 * 1e-6 * T * sqrt(T) / (T + 110.3);
//		lambda = -2 / 3 * mu;
//		k = cp * mu / Pr;
//
//		M[0] = { 0, 0, 0, 0 };
//		M[1][0] = 0;
//		M[1][1] = (lambda + 2 * mu) * nx * nx + mu * ny * ny;
//		M[1][2] = (lambda + mu) * nx * ny;
//		M[1][3] = 0;
//		M[2][0] = 0;
//		M[2][1] = (mu + lambda) * nx * ny;
//		M[2][2] = mu * nx * nx + (lambda + 2 * mu) * ny * ny;
//		M[2][3] = 0;
//		M[3][0] = 0;
//		M[3][1] = u * (lambda + 2 * mu) * nx * nx + (v * mu + v * lambda) * nx * ny + u * mu * ny * ny;
//		M[3][2] = v * mu * nx * nx + (u * lambda + u * mu) * nx * ny + v * (lambda + 2 * mu) * ny * ny;
//		M[3][3] = k * (nx * nx + ny * ny);
//
//		N = {
//			{ {1, 0, 0, 0},
//			{-u / rho, 1 / rho, 0, 0},
//			{-v / rho, 0, 1 / rho, 0},
//			{ (u * u + v * v) / (2 * cv * rho) - T / rho, -u / (cv * rho), -v / (cv * rho), 1 / (cv * rho)}
//			}
//		};
//
//		dx = grid.Center(i, Ny - 1).x - grid.Center(i, Ny - 2).x;
//		dy = grid.Center(i, Ny - 1).y - grid.Center(i, Ny - 2).y;
//
//		dn = sqrt(dx * dx + dy * dy);
//		dv = constoViscPrim(Uii) - constoViscPrim(Ui);
//
//		lp = 0.5 * (Si.uprime + fabs(Si.uprime));
//		lcp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - lp);
//		ltp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));
//
//		V1_Plus = { lcp * Si.k,
//				(Si.u * lcp + Si.a * nx * ltp) * Si.k,
//				(Si.v * lcp + Si.a * ny * ltp) * Si.k,
//				(Si.h0 * lcp + Si.a * Si.uprime * ltp) * Si.k };
//
//		V2_Plus = { ltp / Si.a,
//			Si.u * ltp / Si.a + nx * lcp,
//			Si.v * ltp / Si.a + ny * lcp,
//			Si.h0 * ltp / Si.a + Si.uprime * lcp };
//
//		m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
//		n = { -Si.uprime, nx, ny, 0 };
//
//		X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(4);
//
//		m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
//		n = { -Sii.uprime, nx, ny, 0 };
//
//		lm = 0.5 * (Sii.uprime - fabs(Sii.uprime));
//		lcm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - lm);
//		ltm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));
//
//		V1_Minus = { lcm * Sii.k,
//				(Sii.u * lcm + Sii.a * nx * ltm) * Sii.k,
//				(Sii.v * lcm + Sii.a * ny * ltm) * Sii.k,
//				(Sii.h0 * lcm + Sii.a * Sii.uprime * ltm) * Sii.k };
//
//		V2_Minus = { ltm / Sii.a,
//			Sii.u * ltm / Sii.a + nx * lcm,
//			Sii.v * ltm / Sii.a + ny * lcm,
//			Sii.h0 * ltm / Sii.a + Sii.uprime * lcm };
//
//		Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(4);
//
//		j_viscous_Jacobians[i][Ny] = M / (-dn) * (Etv - identity(4)) * N;
//		j_plus_inviscid_Jacobians[i][Ny] = X;
//		j_minus_inviscid_Jacobians[i][Ny] = Y;
//		j_Fluxes[i][Ny] = X * Ui + Y * Uii - M * dv / dn;
//
//	}
//
//	// Inner i-faces  
//	for (int i = 1; i < Nx; ++i) {
//		for (int j = 0; j < Ny; ++j) {
//
//			nx = grid.iNorms(i, j).x;
//			ny = grid.iNorms(i, j).y;
//
//			Ui = U[i - 1][j];
//			Uii = U[i][j];
//
//			pi = computePressure(Ui);
//			pii = computePressure(Uii);
//			dp = fabs(pii - pi) / min(pii, pi);
//			weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
//
//
//			Up = weight * Ui + (1 - weight) * Uii;
//			Um = (1 - weight) * Ui + weight * Uii;
//
//			Si = compute_viscous_state(Up, nx, ny);
//			Sii = compute_viscous_state(Um, nx, ny);
//
//			W = 0.5 * (Ui + Uii);
//			rho = W[0];
//			u = W[1] / W[0];
//			v = W[2] / W[0];
//			T = computeTemperature(W);
//			mu = 1.458 * 1e-6 * T * sqrt(T) / (T + 110.3);
//			lambda = -2 / 3 * mu;
//			k = cp * mu / Pr;
//
//			M[0] = { 0, 0, 0, 0 };
//			M[1][0] = 0;
//			M[1][1] = (lambda + 2 * mu) * nx * nx + mu * ny * ny;
//			M[1][2] = (lambda + mu) * nx * ny;
//			M[1][3] = 0;
//			M[2][0] = 0;
//			M[2][1] = (mu + lambda) * nx * ny;
//			M[2][2] = mu * nx * nx + (lambda + 2 * mu) * ny * ny;
//			M[2][3] = 0;
//			M[3][0] = 0;
//			M[3][1] = u * (lambda + 2 * mu) * nx * nx + (v * mu + v * lambda) * nx * ny + u * mu * ny * ny;
//			M[3][2] = v * mu * nx * nx + (u * lambda + u * mu) * nx * ny + v * (lambda + 2 * mu) * ny * ny;
//			M[3][3] = k * (nx * nx + ny * ny);
//
//			N = {
//				{ {1, 0, 0, 0},
//				{-u / rho, 1 / rho, 0, 0},
//				{-v / rho, 0, 1 / rho, 0},
//				{ (u * u + v * v) / (2 * cv * rho) - T / rho, -u / (cv * rho), -v / (cv * rho), 1 / (cv * rho)}
//				}
//			};
//
//			dx = grid.Center(i, j).x - grid.Center(i - 1, j).x;
//			dy = grid.Center(i, j).y - grid.Center(i - 1, j).y;
//
//			dn = sqrt(dx * dx + dy * dy);
//			dv = constoViscPrim(Uii) - constoViscPrim(Ui);
//
//			lp = 0.5 * (Si.uprime + fabs(Si.uprime));
//			lcp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - lp);
//			ltp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));
//
//			V1_Plus = { lcp * Si.k,
//				(Si.u * lcp + Si.a * nx * ltp) * Si.k,
//				(Si.v * lcp + Si.a * ny * ltp) * Si.k,
//				(Si.h0 * lcp + Si.a * Si.uprime * ltp) * Si.k };
//
//			V2_Plus = { ltp / Si.a,
//				Si.u * ltp / Si.a + nx * lcp,
//				Si.v * ltp / Si.a + ny * lcp,
//				Si.h0 * ltp / Si.a + Si.uprime * lcp };
//
//			m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
//			n = { -Si.uprime, nx, ny, 0 };
//
//			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(4);
//
//
//			m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
//			n = { -Sii.uprime, nx, ny, 0 };
//
//			lm = 0.5 * (Sii.uprime - fabs(Sii.uprime));
//			lcm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - lm);
//			ltm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));
//
//			V1_Minus = { lcm * Sii.k,
//				(Sii.u * lcm + Sii.a * nx * ltm) * Sii.k,
//				(Sii.v * lcm + Sii.a * ny * ltm) * Sii.k,
//				(Sii.h0 * lcm + Sii.a * Sii.uprime * ltm) * Sii.k };
//
//			V2_Minus = { ltm / Sii.a,
//				Sii.u * ltm / Sii.a + nx * lcm,
//				Sii.v * ltm / Sii.a + ny * lcm,
//				Sii.h0 * ltm / Sii.a + Sii.uprime * lcm };
//
//			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(4);
//
//			i_viscous_Jacobians[i][j] = M / (-dn) * N;
//			i_plus_inviscid_Jacobians[i][j] = X;
//			i_minus_inviscid_Jacobians[i][j] = Y;
//			i_Fluxes[i][j] = X * Ui + Y * Uii - M * dv / dn;
//		}
//	}
//
//	// Inner j-faces 
//	for (int i = 0; i < Nx; ++i) {
//		for (int j = 1; j < Ny; ++j) {
//
//			nx = grid.jNorms(i, j).x;
//			ny = grid.jNorms(i, j).y;
//
//			Ui = U[i][j - 1];
//			Uii = U[i][j];
//
//			pi = computePressure(Ui);
//			pii = computePressure(Uii);
//			dp = fabs(pii - pi) / min(pii, pi);
//			weight = 1 - 0.5 * (1 / ((g * dp) * (g * dp) + 1));
//
//			Up = weight * Ui + (1 - weight) * Uii;
//			Um = (1 - weight) * Ui + weight * Uii;
//
//			Si = compute_viscous_state(Up, nx, ny);
//			Sii = compute_viscous_state(Um, nx, ny);
//
//			W = 0.5 * (Ui + Uii);
//			rho = W[0];
//			u = W[1] / W[0];
//			v = W[2] / W[0];
//			T = computeTemperature(W);
//			mu = 1.458 * 1e-6 * T * sqrt(T) / (T + 110.3);
//			lambda = -2 / 3 * mu;
//			k = cp * mu / Pr;
//
//			M[0] = { 0, 0, 0, 0 };
//			M[1][0] = 0;
//			M[1][1] = (lambda + 2 * mu) * nx * nx + mu * ny * ny;
//			M[1][2] = (lambda + mu) * nx * ny;
//			M[1][3] = 0;
//			M[2][0] = 0;
//			M[2][1] = (mu + lambda) * nx * ny;
//			M[2][2] = mu * nx * nx + (lambda + 2 * mu) * ny * ny;
//			M[2][3] = 0;
//			M[3][0] = 0;
//			M[3][1] = u * (lambda + 2 * mu) * nx * nx + (v * mu + v * lambda) * nx * ny + u * mu * ny * ny;
//			M[3][2] = v * mu * nx * nx + (u * lambda + u * mu) * nx * ny + v * (lambda + 2 * mu) * ny * ny;
//			M[3][3] = k * (nx * nx + ny * ny);
//
//
//			N = {
//				{ {1, 0, 0, 0},
//				{-u / rho, 1 / rho, 0, 0},
//				{-v / rho, 0, 1 / rho, 0},
//				{ (u * u + v * v) / (2 * cv * rho) - T / rho, -u / (cv * rho), -v / (cv * rho), 1 / (cv * rho)}
//				}
//			};
//
//			dx = grid.Center(i, j).x - grid.Center(i, j - 1).x;
//			dy = grid.Center(i, j).y - grid.Center(i, j - 1).y;
//
//			dn = sqrt(dx * dx + dy * dy);
//			dv = constoViscPrim(Uii) - constoViscPrim(Ui);
//
//			lp = 0.5 * (Si.uprime + fabs(Si.uprime));
//			lcp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) + 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)) - lp);
//			ltp = 0.5 * (0.5 * (Si.uprime + Si.a + fabs(Si.uprime + Si.a)) - 0.5 * (Si.uprime - Si.a + fabs(Si.uprime - Si.a)));
//
//			V1_Plus = { lcp * Si.k,
//				(Si.u * lcp + Si.a * nx * ltp) * Si.k,
//				(Si.v * lcp + Si.a * ny * ltp) * Si.k,
//				(Si.h0 * lcp + Si.a * Si.uprime * ltp) * Si.k };
//
//			V2_Plus = { ltp / Si.a,
//				Si.u * ltp / Si.a + nx * lcp,
//				Si.v * ltp / Si.a + ny * lcp,
//				Si.h0 * ltp / Si.a + Si.uprime * lcp };
//
//			m = { Si.pp, -Si.u * pe, -Si.v * pe, pe };
//			n = { -Si.uprime, nx, ny, 0 };
//
//			X = outerProduct(V1_Plus, m) + outerProduct(V2_Plus, n) + lp * identity(4);
//
//			m = { Sii.pp, -Sii.u * pe, -Sii.v * pe, pe };
//			n = { -Sii.uprime, nx, ny, 0 };
//
//			lm = 0.5 * (Sii.uprime - fabs(Sii.uprime));
//			lcm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) + 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)) - lm);
//			ltm = 0.5 * (0.5 * (Sii.uprime + Sii.a - fabs(Sii.uprime + Sii.a)) - 0.5 * (Sii.uprime - Sii.a - fabs(Sii.uprime - Sii.a)));
//
//			V1_Minus = { lcm * Sii.k,
//				(Sii.u * lcm + Sii.a * nx * ltm) * Sii.k,
//				(Sii.v * lcm + Sii.a * ny * ltm) * Sii.k,
//				(Sii.h0 * lcm + Sii.a * Sii.uprime * ltm) * Sii.k };
//
//			V2_Minus = { ltm / Sii.a,
//				Sii.u * ltm / Sii.a + nx * lcm,
//				Sii.v * ltm / Sii.a + ny * lcm,
//				Sii.h0 * ltm / Sii.a + Sii.uprime * lcm };
//
//			Y = outerProduct(V1_Minus, m) + outerProduct(V2_Minus, n) + lm * identity(4);
//
//			j_viscous_Jacobians[i][j] = M / (-dn) * N;
//			j_plus_inviscid_Jacobians[i][j] = X;
//			j_minus_inviscid_Jacobians[i][j] = Y;
//			j_Fluxes[i][j] = X * Ui + Y * Uii - M * dv / dn;
//		}
//	}
//
//}
//void Solver::solve_left_line_viscous() {
//
//	int i = 0;
//	Matrix Eti(4, Vector(4)), Ebi(4, Vector(4)), Eli(4, Vector(4));
//	static Matrix A(4, Vector(4)); // Relates to U(j+1) 
//	static Matrix B(4, Vector(4)); // Relates to U(j) 
//	static Matrix C(4, Vector(4)); // Relates to U(j-1) 
//	static Vector F(4); // Right hand side  
//
//	// Intermediate matrices
//	static Matrix alpha(4, Vector(4));
//	static Matrix v(Ny, Vector(4));
//	static Tensor g(Ny, Matrix(4, Vector(4)));
//
//	// Grab boundary values
//	Ebi = inviscid_boundary_2D_E(BoundaryType.bottom, U[i][0], grid.jNorms(i, 0));
//	Eti = inviscid_boundary_2D_E(BoundaryType.top, U[i][Ny - 1], grid.jNorms(i, Ny));
//
//	Eli = inviscid_boundary_2D_E(BoundaryType.left, U[i][Ny - 1], grid.iNorms(i, Ny - 1));
//
//	// Top Boundary
//	A = grid.Volume(i, Ny - 1) / dt * identity(4)
//		- (i_minus_inviscid_Jacobians[i][Ny - 1] + Eli * i_plus_inviscid_Jacobians[i][Ny - 1] + i_viscous_Jacobians[i][Ny - 1]) * grid.iArea(i, Ny - 1)		// Left     
//		+ (i_plus_inviscid_Jacobians[i + 1][Ny - 1] - i_viscous_Jacobians[i + 1][Ny - 1]) * grid.iArea(i + 1, Ny - 1)										// Right
//		- (j_minus_inviscid_Jacobians[i][Ny - 1] + j_viscous_Jacobians[i][Ny - 1]) * grid.jArea(i, Ny - 1)		 											// Bottom
//		+ (j_plus_inviscid_Jacobians[i][Ny] + Eti * j_minus_inviscid_Jacobians[i][Ny] + j_viscous_Jacobians[i][Ny]) * grid.jArea(i, Ny);					// Top  
//
//	C = (j_plus_inviscid_Jacobians[i][Ny - 1] - j_viscous_Jacobians[i][Ny - 1]) * (-grid.jArea(i, Ny - 1));
//
//	F = i_Fluxes[i][Ny - 1] * (grid.iArea(i, Ny - 1))
//		+ (i_Fluxes[i + 1][Ny - 1] + (i_minus_inviscid_Jacobians[i + 1][Ny - 1] + i_viscous_Jacobians[i + 1][Ny - 1]) * dU_old[i + 1][Ny - 1]) * (-grid.iArea(i + 1, Ny - 1))
//		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
//		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny));
//
//	alpha = A;
//	v[Ny - 1] = F / alpha;
//	g[Ny - 1] = C / alpha;
//
//	// Middle Boundry
//	for (int j = Ny - 2; j > 0; --j) {
//
//		Eli = inviscid_boundary_2D_E(BoundaryType.left, U[i][j], grid.iNorms(i, j));
//
//		B = (j_minus_inviscid_Jacobians[i][j + 1] + j_viscous_Jacobians[i][j + 1]) * (grid.jArea(i, j + 1));
//
//		A = grid.Volume(i, j) / dt * identity(4)
//			- (i_minus_inviscid_Jacobians[i][j] + Eli * i_plus_inviscid_Jacobians[i][j] + i_viscous_Jacobians[i][j]) * grid.iArea(i, j)			// Left  
//			+ (i_plus_inviscid_Jacobians[i + 1][j] - i_viscous_Jacobians[i + 1][j]) * grid.iArea(i + 1, j)										// Right
//			- (j_minus_inviscid_Jacobians[i][j] + j_viscous_Jacobians[i][j]) * grid.jArea(i, j)													// Bottom
//			+ (j_plus_inviscid_Jacobians[i][j + 1] - j_viscous_Jacobians[i][j + 1]) * grid.jArea(i, j + 1);										// Top
//
//		C = (j_plus_inviscid_Jacobians[i][j] - j_viscous_Jacobians[i][j]) * (-grid.jArea(i, j));
//
//		F = i_Fluxes[i][j] * (grid.iArea(i, j))
//			+ (i_Fluxes[i + 1][j] + (i_minus_inviscid_Jacobians[i + 1][j] + i_viscous_Jacobians[i + 1][j]) * dU_old[i + 1][j]) * (-grid.iArea(i + 1, j))
//			+ j_Fluxes[i][j] * (grid.jArea(i, j))
//			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1));
//
//
//
//		alpha = A - B * g[j + 1];
//		g[j] = C / alpha;
//		v[j] = (F - B * v[j + 1]) / alpha;
//	}
//
//	//  Bottom boundary
//
//	Eli = inviscid_boundary_2D_E(BoundaryType.left, U[i][0], grid.iNorms(i, 0));
//
//	B = (j_minus_inviscid_Jacobians[i][1] + j_viscous_Jacobians[i][1]) * (grid.jArea(i, 1));
//
//	A = grid.Volume(i, 0) / dt * identity(4)
//		- (i_minus_inviscid_Jacobians[i][0] + Eli * i_plus_inviscid_Jacobians[i][0] + i_viscous_Jacobians[i][0]) * grid.iArea(i, 0)		// Left
//		+ (i_plus_inviscid_Jacobians[i + 1][0] - i_viscous_Jacobians[i + 1][0]) * grid.iArea(i + 1, 0)									// Right 
//		- (Ebi * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0] + j_viscous_Jacobians[i][0]) * grid.jArea(i, 0)		// Bottom
//		+ (j_plus_inviscid_Jacobians[i][1] - j_viscous_Jacobians[i][1]) * grid.jArea(i, 1);												// Top
//
//	F = i_Fluxes[i][0] * (grid.iArea(i, 0))
//		+ (i_Fluxes[i + 1][0] + (i_minus_inviscid_Jacobians[i + 1][0] + i_viscous_Jacobians[i + 1][0]) * dU_old[i + 1][0]) * (-grid.iArea(i + 1, 0))
//		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
//		+ j_Fluxes[i][1] * (-grid.jArea(i, 1));
//
//
//	alpha = A - B * g[1];
//	v[0] = (F - B * v[1]) / alpha;
//
//	// Calculate dU
//	dU_new[i][0] = v[0];
//	for (int j = 1; j < Ny; ++j) {
//		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];
//	}
//
//}
//void Solver::solve_middle_line_viscous(const int i) {
//
//	Matrix Eti(4, Vector(4)), Ebi(4, Vector(4));
//
//	static Matrix A(4, Vector(4)); // Relates to U(j+1) 
//	static Matrix B(4, Vector(4)); // Relates to U(j) 
//	static Matrix C(4, Vector(4)); // Relates to U(j-1) 
//	static Vector F(4); // Right hand side  
//
//	// Intermediate matrices
//	static Matrix alpha(4, Vector(4));
//	static Matrix v(Ny, Vector(4));
//	static Tensor g(Ny, Matrix(4, Vector(4)));
//
//	// Grab boundary values
//	Ebi = inviscid_boundary_2D_E(BoundaryType.bottom, U[i][0], grid.jNorms(i, 0));
//	Eti = inviscid_boundary_2D_E(BoundaryType.top, U[i][Ny - 1], grid.jNorms(i, Ny));
//
//
//	// Top Boundary
//	A = grid.Volume(i, Ny - 1) / dt * identity(4)
//		- (i_minus_inviscid_Jacobians[i][Ny - 1] + i_viscous_Jacobians[i][Ny - 1]) * grid.iArea(i, Ny - 1)															// Left
//		+ (i_plus_inviscid_Jacobians[i + 1][Ny - 1] - i_viscous_Jacobians[i + 1][Ny - 1]) * grid.iArea(i + 1, Ny - 1)												// Right
//		- (j_minus_inviscid_Jacobians[i][Ny - 1] + j_viscous_Jacobians[i][Ny - 1]) * grid.jArea(i, Ny - 1)															// Bottom
//		+ (j_plus_inviscid_Jacobians[i][Ny] + Eti * j_minus_inviscid_Jacobians[i][Ny] + j_viscous_Jacobians[i][Ny]) * grid.jArea(i, Ny);		// Right 
//
//	C = (j_plus_inviscid_Jacobians[i][Ny - 1] - j_viscous_Jacobians[i][Ny - 1]) * (-grid.jArea(i, Ny - 1));
//
//	F = (i_Fluxes[i][Ny - 1] + (i_plus_inviscid_Jacobians[i][Ny - 1] - i_viscous_Jacobians[i][Ny - 1]) * dU_old[i - 1][Ny - 1]) * (grid.iArea(i, Ny - 1))
//		+ (i_Fluxes[i + 1][Ny - 1] + (i_minus_inviscid_Jacobians[i + 1][Ny - 1] + i_viscous_Jacobians[i + 1][Ny - 1]) * dU_old[i + 1][Ny - 1]) * (-grid.iArea(i + 1, Ny - 1))
//		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
//		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny));
//
//	alpha = A;
//	v[Ny - 1] = F / alpha;
//	g[Ny - 1] = C / alpha;
//
//	// Middle Boundry
//	for (int j = Ny - 2; j > 0; --j) {
//
//		B = (j_minus_inviscid_Jacobians[i][j + 1] + j_viscous_Jacobians[i][j + 1]) * (grid.jArea(i, j + 1));
//
//		A = grid.Volume(i, j) / dt * identity(4)
//			- (i_minus_inviscid_Jacobians[i][j] + i_viscous_Jacobians[i][j]) * grid.iArea(i, j)
//			+ (i_plus_inviscid_Jacobians[i + 1][j] - i_viscous_Jacobians[i + 1][j]) * grid.iArea(i + 1, j)
//			- (j_minus_inviscid_Jacobians[i][j] + j_viscous_Jacobians[i][j]) * grid.jArea(i, j)
//			+ (j_plus_inviscid_Jacobians[i][j + 1] - j_viscous_Jacobians[i][j + 1]) * grid.jArea(i, j + 1);
//
//		C = (j_plus_inviscid_Jacobians[i][j] - j_viscous_Jacobians[i][j]) * (-grid.jArea(i, j));
//
//		F = (i_Fluxes[i][j] + (i_plus_inviscid_Jacobians[i][j] - i_viscous_Jacobians[i][j]) * dU_old[i - 1][j]) * (grid.iArea(i, j))
//			+ (i_Fluxes[i + 1][j] + (i_minus_inviscid_Jacobians[i + 1][j] + i_viscous_Jacobians[i + 1][j]) * dU_old[i + 1][j]) * (-grid.iArea(i + 1, j))
//			+ j_Fluxes[i][j] * (grid.jArea(i, j))
//			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1));
//
//
//
//		alpha = A - B * g[j + 1];
//		g[j] = C / alpha;
//		v[j] = (F - B * v[j + 1]) / alpha;
//	}
//
//	//  Bottom boundary
//	B = (j_minus_inviscid_Jacobians[i][1] + j_viscous_Jacobians[i][1]) * (grid.jArea(i, 1));
//
//	A = grid.Volume(i, 0) / dt * identity(4)
//		- (i_minus_inviscid_Jacobians[i][0] + i_viscous_Jacobians[i][0]) * grid.iArea(i, 0)
//		+ (i_plus_inviscid_Jacobians[i + 1][0] - i_viscous_Jacobians[i + 1][0]) * grid.iArea(i + 1, 0)
//		- (Ebi * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0] + j_viscous_Jacobians[i][0]) * grid.jArea(i, 0)
//		+ (j_plus_inviscid_Jacobians[i][1] - j_viscous_Jacobians[i][1]) * grid.jArea(i, 1);
//
//	F = (i_Fluxes[i][0] + (i_plus_inviscid_Jacobians[i][0] - i_viscous_Jacobians[i][0]) * dU_old[i - 1][0]) * (grid.iArea(i, 0))
//		+ (i_Fluxes[i + 1][0] + (i_minus_inviscid_Jacobians[i + 1][0] + i_viscous_Jacobians[i + 1][0]) * dU_old[i + 1][0]) * (-grid.iArea(i + 1, 0))
//		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
//		+ j_Fluxes[i][1] * (-grid.jArea(i, 1));
//
//
//	alpha = A - B * g[1];
//	v[0] = (F - B * v[1]) / alpha;
//
//	// Calculate dU
//	dU_new[i][0] = v[0];
//
//	for (int j = 1; j < Ny; ++j) {
//		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];
//	}
//}
//void Solver::solve_right_line_viscous() {
//
//	int i = Nx - 1;
//	Matrix Eti(4, Vector(4)), Ebi(4, Vector(4)), Eri(4, Vector(4));
//
//	static Matrix A(4, Vector(4)); // Relates to U(j+1) 
//	static Matrix B(4, Vector(4)); // Relates to U(j) 
//	static Matrix C(4, Vector(4)); // Relates to U(j-1) 
//	static Vector F(4); // Right hand side  
//
//	// Intermediate matrices
//	static Matrix alpha(4, Vector(4));
//	static Matrix v(Ny, Vector(4));
//	static Tensor g(Ny, Matrix(4, Vector(4)));
//
//	// Grab boundary values
//	Ebi = inviscid_boundary_2D_E(BoundaryType.bottom, U[i][0], grid.jNorms(i, 0));
//	Eti = inviscid_boundary_2D_E(BoundaryType.top, U[i][Ny - 1], grid.jNorms(i, Ny));
//
//
//	Eri = inviscid_boundary_2D_E(BoundaryType.right, U[i][Ny - 1], grid.iNorms(i + 1, Ny - 1));
//
//	// Top Boundary
//	A = grid.Volume(i, Ny - 1) / dt * identity(4)
//		- (i_minus_inviscid_Jacobians[i][Ny - 1] + i_viscous_Jacobians[i][Ny - 1]) * grid.iArea(i, Ny - 1)																// Left
//		+ (i_plus_inviscid_Jacobians[i + 1][Ny - 1] + Eri * i_minus_inviscid_Jacobians[i + 1][Ny - 1] + i_viscous_Jacobians[i + 1][Ny - 1]) * grid.iArea(i + 1, Ny - 1)	// Right
//		- (j_minus_inviscid_Jacobians[i][Ny - 1] + j_viscous_Jacobians[i][Ny - 1]) * grid.jArea(i, Ny - 1)																// Bottom
//		+ (j_plus_inviscid_Jacobians[i][Ny] + Eti * j_minus_inviscid_Jacobians[i][Ny] + j_viscous_Jacobians[i][Ny]) * grid.jArea(i, Ny);		 						// Top 
//
//	C = (j_plus_inviscid_Jacobians[i][Ny - 1] - j_viscous_Jacobians[i][Ny - 1]) * (-grid.jArea(i, Ny - 1));
//
//	F = (i_Fluxes[i][Ny - 1] + (i_plus_inviscid_Jacobians[i][Ny - 1] - i_viscous_Jacobians[i][Ny - 1]) * dU_old[i - 1][Ny - 1]) * (grid.iArea(i, Ny - 1))
//		+ i_Fluxes[i + 1][Ny - 1] * (-grid.iArea(i + 1, Ny - 1))
//		+ j_Fluxes[i][Ny - 1] * (grid.jArea(i, Ny - 1))
//		+ j_Fluxes[i][Ny] * (-grid.jArea(i, Ny));
//
//	alpha = A;
//	v[Ny - 1] = F / alpha;
//	g[Ny - 1] = C / alpha;
//
//	// Middle Boundry
//	for (int j = Ny - 2; j > 0; --j) {
//
//		Eri = inviscid_boundary_2D_E(BoundaryType.right, U[i][j], grid.iNorms(i + 1, j));
//
//		B = (j_minus_inviscid_Jacobians[i][j + 1] + j_viscous_Jacobians[i][j + 1]) * (grid.jArea(i, j + 1));
//
//		A = grid.Volume(i, j) / dt * identity(4)
//			- (i_minus_inviscid_Jacobians[i][j] + i_viscous_Jacobians[i][j]) * grid.iArea(i, j)															// Left
//			+ (i_plus_inviscid_Jacobians[i + 1][j] + Eri * i_minus_inviscid_Jacobians[i + 1][j] + i_viscous_Jacobians[i + 1][j]) * grid.iArea(i + 1, j)	// Right
//			- (j_minus_inviscid_Jacobians[i][j] + j_viscous_Jacobians[i][j]) * grid.jArea(i, j)															// Bottom
//			+ (j_plus_inviscid_Jacobians[i][j + 1] - j_viscous_Jacobians[i][j + 1]) * grid.jArea(i, j + 1);												// Top
//
//		C = (j_plus_inviscid_Jacobians[i][j] - j_viscous_Jacobians[i][j]) * (-grid.jArea(i, j));
//
//		F = (i_Fluxes[i][j] + (i_plus_inviscid_Jacobians[i][j] - i_viscous_Jacobians[i][j]) * dU_old[i - 1][j]) * (grid.iArea(i, j))
//			+ i_Fluxes[i + 1][j] * (-grid.iArea(i + 1, j))
//			+ j_Fluxes[i][j] * (grid.jArea(i, j))
//			+ j_Fluxes[i][j + 1] * (-grid.jArea(i, j + 1));
//
//		alpha = A - B * g[j + 1];
//		g[j] = C / alpha;
//		v[j] = (F - B * v[j + 1]) / alpha;
//	}
//	//  Bottom boundary
//	Eri = inviscid_boundary_2D_E(BoundaryType.right, U[i][0], grid.iNorms(i + 1, 0));
//
//	B = (j_minus_inviscid_Jacobians[i][1] + j_viscous_Jacobians[i][1]) * (grid.jArea(i, 1));
//
//	A = grid.Volume(i, 0) / dt * identity(4)
//		- (i_minus_inviscid_Jacobians[i][0] + i_viscous_Jacobians[i][0]) * grid.iArea(i, 0)															// Left
//		+ (i_plus_inviscid_Jacobians[i + 1][0] + Eri * i_minus_inviscid_Jacobians[i + 1][0] + i_viscous_Jacobians[i + 1][0]) * grid.iArea(i + 1, 0)	// Right 
//		- (Ebi * j_plus_inviscid_Jacobians[i][0] + j_minus_inviscid_Jacobians[i][0] + j_viscous_Jacobians[i][0]) * grid.jArea(i, 0)					// Bottom 
//		+ (j_plus_inviscid_Jacobians[i][1] - j_viscous_Jacobians[i][1]) * grid.jArea(i, 1);															// Top 
//
//	F = (i_Fluxes[i][0] + (i_plus_inviscid_Jacobians[i][0] - i_viscous_Jacobians[i][0]) * dU_old[i - 1][0]) * (grid.iArea(i, 0))
//		+ i_Fluxes[i + 1][0] * (-grid.iArea(i + 1, 0))
//		+ j_Fluxes[i][0] * (grid.jArea(i, 0))
//		+ j_Fluxes[i][1] * (-grid.jArea(i, 1));
//
//
//	alpha = A - B * g[1];
//	v[0] = (F - B * v[1]) / alpha;
//
//	// Calculate dU
//	dU_new[i][0] = v[0];
//
//	for (int j = 1; j < Ny; ++j) {
//		dU_new[i][j] = v[j] - g[j] * dU_new[i][j - 1];
//	}
//
//}
//void Solver::viscous_calculations() {
//
//	Vector tau(Nx), q(Nx), y_plus(Nx), Re_cell(Nx), velocity_profile(Ny), temperature_profile(Ny);
//	double rho, u, T, mu, k, ut, dy;
//
//	for (int i = 0; i < Nx; ++i) {
//
//		rho = U[i][0][0];
//		u = U[i][0][1] / U[i][0][0];
//		T = computeTemperature(U[i][0]);
//		mu = 1.458 * 1e-6 * T * sqrt(T) / (T + 110.3);
//		k = cp * mu / Pr;
//		dy = grid.Center(i, 0).y;
//
//
//		tau[i] = mu * u / dy;
//		q[i] = -k * (T - Tw) / dy;
//		y_plus[i] = rho * sqrt(tau[i] / rho) * dy / mu;
//		Re_cell[i] = rho * u * grid.jArea(i, 0) / mu;
//
//	}
//
//
//	for (int j = 0; j < Ny; ++j) {
//		u = U[0][j][1] / U[0][j][0];
//		T = computeTemperature(U[0][j]);
//
//		velocity_profile[j] = u;
//		temperature_profile[j] = T;
//	}
//
//	string filename = "Viscous_gradient_inclusive_terms.csv";
//	ofstream file(filename);
//
//	file << "x_points, tau, q, y_plus, Re_cell";
//	for (int i = 0; i < Nx; ++i) {
//		file << grid.Center(i, 0).x << ", " << tau[i] << ", " << q[i] << ", " << y_plus[i] << ", " << Re_cell[i] << endl;
//	}
//
//	file.close();
//
//
//	filename = "Velocity/Temperature_profile.csv";
//	ofstream file2(filename);
//	file2 << "y_points, velocity_profile, temperature_profile" << endl;
//	for (int j = 0; j < Ny; ++j) {
//		file2 << grid.Center(0, j).y << ", " << velocity_profile[j] << ", " << temperature_profile[j] << endl;
//	}
//	file2.close();
//
//
//}
