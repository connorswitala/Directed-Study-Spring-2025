#pragma once

#include "pch.h"
#include "framework.h"
#include <iostream>
#include "LinearAlgebra.h"
#include <fstream>

using namespace std;

struct Point {
    double x;
    double y;

    Point(double x_val = 0.0, double y_val = 0.0) : x(x_val), y(y_val) {};
};



class Grid {
public:
    virtual double Volume(int i, int j) const = 0;
    virtual double LArea(int i, int j) const = 0;
    virtual double RArea(int i, int j) const = 0;
    virtual double BArea(int i, int j) const = 0;
    virtual double TArea(int i, int j) const = 0;
    virtual Point LNorms(int i, int j) const = 0;
    virtual Point RNorms(int i, int j) const = 0;
    virtual Point BNorms(int i, int j) const = 0;
    virtual Point TNorms(int i, int j) const = 0;
    virtual Point Center(int i, int j) const = 0;
    virtual Point Vertex(int i, int j) const = 0;


    virtual ~Grid() = default;
};


class RampGrid : public Grid {
private:
    int Nx, Ny;
    double L1, L2, L3, inlet_height, ramp_angle;

    vector<vector<Point>> vertices;
    vector<vector<Point>> cellCenters;
    vector<vector<Point>> LNormals;
    vector<vector<Point>> RNormals;
    vector<vector<Point>> BNormals;
    vector<vector<Point>> TNormals;

    Matrix LAreas, RAreas, BAreas, TAreas, cellVolumes;

public:
    RampGrid(int Nx, int Ny, double L1, double L2, double L3, double inlet_height, double ramp_angle);

    double Volume(int i, int j) const override;
    double LArea(int i, int j) const override;
    double RArea(int i, int j) const override;
    double BArea(int i, int j) const override;
    double TArea(int i, int j) const override;
    Point LNorms(int i, int j) const override;
    Point RNorms(int i, int j) const override;
    Point BNorms(int i, int j) const override;
    Point TNorms(int i, int j) const override;
    Point Center(int i, int j) const override;
    Point Vertex(int i, int j) const override;

};

class SquareGrid : public Grid {
private:
    int Lx, Nx, Ly, Ny;
    vector<vector<Point>> vertices;
    vector<vector<Point>> cellCenters;
    vector<vector<Point>> faceAreas;
    vector<vector<Point>> LNormals;
    vector<vector<Point>> RNormals;
    vector<vector<Point>> BNormals;
    vector<vector<Point>> TNormals;
    Matrix LAreas;
    Matrix RAreas;
    Matrix BAreas;
    Matrix TAreas;
    Matrix cellVolumes;

public:

    SquareGrid(int Lx, int Nx, int Ly, int Ny);

    double Volume(int i, int j) const override;
    double LArea(int i, int j) const override;
    double RArea(int i, int j) const override;
    double BArea(int i, int j) const override;
    double TArea(int i, int j) const override;
    Point LNorms(int i, int j) const override;
    Point RNorms(int i, int j) const override;
    Point BNorms(int i, int j) const override;
    Point TNorms(int i, int j) const override;
    Point Center(int i, int j) const override;
    Point Vertex(int i, int j) const override;




};