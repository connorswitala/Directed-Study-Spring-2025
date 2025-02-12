
#include "pch.h"
#include "framework.h"
#include "GridGenerator.h"

/////////////////////////////////////////////////
///////////// Ramp Grid functions ///////////////
/////////////////////////////////////////////////

RampGrid::RampGrid(int Nx, int Ny, double L1, double L2, double L3, double inlet_height, double ramp_angle)
	: Nx(Nx), Ny(Ny), L1(L1), L2(L2), L3(L3), inlet_height(inlet_height), ramp_angle(ramp_angle),
	vertices(Nx + 1, vector<Point>(Ny + 1)), cellCenters(Nx, vector<Point>(Ny)),
	cellVolumes(Nx, vector<double>(Ny)), LAreas(Nx, Vector(Ny, 0.0)),
	RAreas(Nx, Vector(Ny, 0.0)), BAreas(Nx, Vector(Ny, 0.0)), TAreas(Nx, Vector(Ny, 0.0)),
	LNormals(Nx, vector<Point>(Ny)), RNormals(Nx, vector<Point>(Ny)), BNormals(Nx, vector<Point>(Ny)),
	TNormals(Nx, vector<Point>(Ny)) {


	int i, j;


	// Important constants
	double deg_to_rads = 3.141592653 / 180.0;
	double ang_rads = ramp_angle * deg_to_rads;
	double min_height = L2 * tan(ang_rads) / 0.9;

	// Checks that inlet height is tall enough for ramp section
	if (inlet_height < min_height) {
		cout << "Ramp angle is too large for inlet height, please input an inlet height larger than " << min_height << endl;
		cin >> inlet_height;
	}

	// Define lengths 
	double L_total = L1 + L2 + L3;
	double dx = L_total / (Nx + 1);
	double dy = inlet_height / (Ny + 1);
	double L, dy_ramp;

	// Create x vertices
	for (i = 0; i <= Nx; ++i) {
		for (j = 0; j <= Ny; ++j) {
			vertices[i][j].x = i * dx;
		}
	}

	// Snap x-vertices to important boundary points.
	for (i = 0; i <= Nx; ++i) {
		for (j = 0; j <= Ny; ++j) {
			if (vertices[i][j].x > L1 && vertices[i][j].x < L1 + dx) {
				vertices[i][j].x = L1;
			}
			else if (vertices[i][j].x > L1 + L2 && vertices[i][j].x < L1 + L2 + dx) {
				vertices[i][j].x = L1 + L2;
			}
			else if (i == Nx) {
				vertices[i][j].x = L_total;
			}
		}
	}


	// Create grid.
	for (i = 0; i <= Nx; ++i) {
		L = vertices[i][0].x;

		// y vertices for section 1
		if (L <= L1) {
			dy = inlet_height / Ny;
			for (j = 0; j <= Ny; ++j) {
				vertices[i][j].y = j * dy;
			}
		}

		// y vertices for section 2
		else if ((L > L1) && (L <= (L1 + L2))) {

			// Changes dy based on section
			if (vertices[i - 1][0].x == L1) {
				dy_ramp = (vertices[i][0].x - vertices[i - 1][0].x) * tan(ang_rads);
			}
			else if ((L + dx > L1 + L2) && (L < L1 + L2 + dx)) {
				dy_ramp = (vertices[i][0].x - vertices[i - 1][0].x) * tan(ang_rads);
			}
			else {
				dy_ramp = dx * tan(ang_rads);
			}

			vertices[i][0].y = vertices[i - 1][0].y + dy_ramp;
			dy = (inlet_height - vertices[i][0].y) / Ny;

			for (j = 1; j <= Ny; ++j) {
				vertices[i][j].y = vertices[i][0].y + j * dy;
			}
		}

		// y vertives for section 3
		else {
			vertices[i][0].y = L2 * tan(ang_rads);
			dy = (inlet_height - L2 * tan(ang_rads)) / Ny;

			for (j = 1; j <= Ny; ++j) {
				vertices[i][j].y = L2 * tan(ang_rads) + j * dy;
			}
		}
	}

	// Edge vectors
	vector<vector<Point>> DA(Nx, vector<Point>(Ny)), AB(Nx, vector<Point>(Ny)), BC(Nx, vector<Point>(Ny)), CD(Nx, vector<Point>(Ny));

	// Calculate i-face normals
	for (i = 0; i < Nx; ++i) {
		for (j = 0; j < Ny; ++j) {
			DA[i][j] = { vertices[i][j].x - vertices[i][j + 1].x, vertices[i][j].y - vertices[i][j + 1].y };
			AB[i][j] = { vertices[i + 1][j].x - vertices[i][j].x, vertices[i + 1][j].y - vertices[i][j].y };
			BC[i][j] = { vertices[i + 1][j + 1].x - vertices[i + 1][j].x, vertices[i + 1][j + 1].y - vertices[i + 1][j].y };
			CD[i][j] = { vertices[i][j + 1].x - vertices[i + 1][j + 1].x, vertices[i][j + 1].y - vertices[i + 1][j + 1].y };


			LAreas[i][j] = sqrt(DA[i][j].x * DA[i][j].x + DA[i][j].y * DA[i][j].y);
			BAreas[i][j] = sqrt(AB[i][j].x * AB[i][j].x + AB[i][j].y * AB[i][j].y);
			RAreas[i][j] = sqrt(BC[i][j].x * BC[i][j].x + BC[i][j].y * BC[i][j].y);
			TAreas[i][j] = sqrt(CD[i][j].x * CD[i][j].x + CD[i][j].y * CD[i][j].y);

			LNormals[i][j].x = DA[i][j].y / fabs(LAreas[i][j]);
			LNormals[i][j].y = DA[i][j].x / fabs(LAreas[i][j]);

			BNormals[i][j].x = AB[i][j].y / fabs(BAreas[i][j]);
			BNormals[i][j].y = -AB[i][j].x / fabs(BAreas[i][j]);

			RNormals[i][j].x = BC[i][j].y / fabs(RAreas[i][j]);
			RNormals[i][j].y = BC[i][j].x / fabs(RAreas[i][j]);

			TNormals[i][j].x = CD[i][j].y / fabs(TAreas[i][j]);
			TNormals[i][j].y = -CD[i][j].x / fabs(TAreas[i][j]);

			cellCenters[i][j] = { (vertices[i][j].x + vertices[i + 1][j].x + vertices[i + 1][j + 1].x + vertices[i][j + 1].x) / 4, (vertices[i][j].y + vertices[i + 1][j].y + vertices[i + 1][j + 1].y + vertices[i][j + 1].y) / 4 };
			cellVolumes[i][j] = 0.5 * fabs(DA[i][j].x * AB[i][j].y - DA[i][j].y * AB[i][j].x) + 0.5 * fabs(BC[i][j].x * CD[i][j].y - BC[i][j].y * CD[i][j].x);


		}
	}
}



Point RampGrid::Center(int i, int j) const {
	return cellCenters[i][j];
}

Point RampGrid::Vertex(int i, int j) const {
	return vertices[i][j];
}

double RampGrid::Volume(int i, int j) const {
	return cellVolumes[i][j];
}

double RampGrid::LArea(int i, int j) const {
	return LAreas[i][j];
}

double RampGrid::RArea(int i, int j) const {
	return RAreas[i][j];
}

double RampGrid::BArea(int i, int j) const {
	return BAreas[i][j];
}

double RampGrid::TArea(int i, int j) const {
	return TAreas[i][j];
}

Point RampGrid::LNorms(int i, int j) const {
	return LNormals[i][j];
}
Point RampGrid::RNorms(int i, int j) const {
	return RNormals[i][j];
}
Point RampGrid::BNorms(int i, int j) const {
	return BNormals[i][j];
}
Point RampGrid::TNorms(int i, int j) const {
	return TNormals[i][j];
}

/////////////////////////////////////////////////
///////////// Square Grid functions /////////////
/////////////////////////////////////////////////

SquareGrid::SquareGrid(int Lx, int Nx, int Ly, int Ny)
	: Nx(Nx), Ny(Ny), Lx(Lx), Ly(Ly),
	vertices(Nx + 1, vector<Point>(Ny + 1)), cellCenters(Nx, vector<Point>(Ny)),
	cellVolumes(Nx, vector<double>(Ny)), LAreas(Nx, Vector(Ny, 0.0)),
	RAreas(Nx, Vector(Ny, 0.0)), BAreas(Nx, Vector(Ny, 0.0)), TAreas(Nx, Vector(Ny, 0.0)),
	LNormals(Nx, vector<Point>(Ny)), RNormals(Nx, vector<Point>(Ny)), BNormals(Nx, vector<Point>(Ny)),
	TNormals(Nx, vector<Point>(Ny)) {


	int i, j;


	// Define lengths 
	double dx = Lx / (Nx + 1);
	double dy = Ly / (Ny + 1);

	// Create x vertices
	for (i = 0; i <= Nx; ++i) {
		for (j = 0; j <= Ny; ++j) {
			vertices[i][j].x = i * dx;
			vertices[i][j].y = j * dy;
		}
	}

	// Edge vectors
	vector<vector<Point>> DA(Nx, vector<Point>(Ny)), AB(Nx, vector<Point>(Ny)), BC(Nx, vector<Point>(Ny)), CD(Nx, vector<Point>(Ny));

	// Calculate i-face normals
	for (i = 0; i < Nx; ++i) {
		for (j = 0; j < Ny; ++j) {
			DA[i][j] = { vertices[i][j].x - vertices[i][j + 1].x, vertices[i][j].y - vertices[i][j + 1].y };
			AB[i][j] = { vertices[i + 1][j].x - vertices[i][j].x, vertices[i + 1][j].y - vertices[i][j].y };
			BC[i][j] = { vertices[i + 1][j + 1].x - vertices[i + 1][j].x, vertices[i + 1][j + 1].y - vertices[i + 1][j].y };
			CD[i][j] = { vertices[i][j + 1].x - vertices[i + 1][j + 1].x, vertices[i][j + 1].y - vertices[i + 1][j + 1].y };


			LAreas[i][j] = sqrt(DA[i][j].x * DA[i][j].x + DA[i][j].y * DA[i][j].y);
			BAreas[i][j] = sqrt(AB[i][j].x * AB[i][j].x + AB[i][j].y * AB[i][j].y);
			RAreas[i][j] = sqrt(BC[i][j].x * BC[i][j].x + BC[i][j].y * BC[i][j].y);
			TAreas[i][j] = sqrt(CD[i][j].x * CD[i][j].x + CD[i][j].y * CD[i][j].y);

			LNormals[i][j].x = DA[i][j].y / fabs(LAreas[i][j]);
			LNormals[i][j].y = DA[i][j].x / fabs(LAreas[i][j]);

			BNormals[i][j].x = AB[i][j].y / fabs(BAreas[i][j]);
			BNormals[i][j].y = -AB[i][j].x / fabs(BAreas[i][j]);

			RNormals[i][j].x = BC[i][j].y / fabs(RAreas[i][j]);
			RNormals[i][j].y = BC[i][j].x / fabs(RAreas[i][j]);

			TNormals[i][j].x = CD[i][j].y / fabs(TAreas[i][j]);
			TNormals[i][j].y = -CD[i][j].x / fabs(TAreas[i][j]);

			cellCenters[i][j] = { (vertices[i][j].x + vertices[i + 1][j].x + vertices[i + 1][j + 1].x + vertices[i][j + 1].x) / 4, (vertices[i][j].y + vertices[i + 1][j].y + vertices[i + 1][j + 1].y + vertices[i][j + 1].y) / 4 };
			cellVolumes[i][j] = 0.5 * fabs(DA[i][j].x * AB[i][j].y - DA[i][j].y * AB[i][j].x) + 0.5 * fabs(BC[i][j].x * CD[i][j].y - BC[i][j].y * CD[i][j].x);


		}
	}
}



Point SquareGrid::SquareGrid::Center(int i, int j) const {
	return cellCenters[i][j];
}

Point SquareGrid::SquareGrid::Vertex(int i, int j) const {
	return vertices[i][j];
}

double SquareGrid::SquareGrid::Volume(int i, int j) const {
	return cellVolumes[i][j];
}

double SquareGrid::SquareGrid::LArea(int i, int j) const {
	return LAreas[i][j];
}

double SquareGrid::SquareGrid::RArea(int i, int j) const {
	return RAreas[i][j];
}

double SquareGrid::SquareGrid::BArea(int i, int j) const {
	return BAreas[i][j];
}

double SquareGrid::SquareGrid::TArea(int i, int j) const {
	return TAreas[i][j];
}

Point SquareGrid::SquareGrid::LNorms(int i, int j) const {
	return LNormals[i][j];
}
Point SquareGrid::SquareGrid::RNorms(int i, int j) const {
	return RNormals[i][j];
}
Point SquareGrid::SquareGrid::BNorms(int i, int j) const {
	return BNormals[i][j];
}
Point SquareGrid::SquareGrid::TNorms(int i, int j) const {
	return TNormals[i][j];
}