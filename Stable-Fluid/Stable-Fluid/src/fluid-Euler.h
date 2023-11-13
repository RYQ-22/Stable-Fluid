#ifndef FLUID_EULER_H
#define FLUID_EULER_H

#include "draw.h"
#include "vector3.h"

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include <ctime>

#include <vector>
#include <unordered_map>
#include <queue>
#include <thread>
#include <functional>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#define PI 3.1415926f
#define G_ 9.8f
#define Gamma_ 0.1f

class Grid {
private:
	int N1_;
	int N2_;
	int N3_;
	int N_fluid_ = 0;
	float l_;
	// position of the grid's center
	Vec3<float> x_, y_, z_;
	// velocity at the grid's edge
	Vec3<float> Vx_, Vy_, Vz_;
	// level set
	Vec3<float> phi_;
	// boundary
	Vec3<int> solid_;
	Vec3<int> idx_;
	Vec3<int> Vx_valid_;
	Vec3<int> Vy_valid_;
	Vec3<int> Vz_valid_;
	// matrix
	Vec3<int> Adiag, Aplusi, Aplusj, Aplusk;
	Vec3<float> d, p, precon, q;

public:
	Grid();
	Grid(int N1, int N2, int N3, float l, std::vector<float> phi, std::vector<int> solid);
	float Vx_ave(float i, float j, float k);
	float Vy_ave(float i, float j, float k);
	float Vz_ave(float i, float j, float k);
	float phi_ave(float i, float j, float k);

	friend class Fluid_Euler;
};

class Fluid_Euler {
private:
	float nu_ = 0.1f;
	float l_;
	float rho_ = 1.0f;
	Grid grid_;
public:
	Fluid_Euler(int N1, int N2, int N3, float l, std::vector<float> phi, std::vector<int> solid);

	bool valid(int i, int j, int k);// test if (i, j, k) is liquid

	void update(float dt);

	// update fluid
	void add_force(float dt);// 1. add force
	void advect(float dt);// 2. advect
	void project(float dt);// 3. project
	void solve();// solve: A * p = d;

	// some details
	void compute_phi(float dt);
	void extrapolate();
	void constrain();

	// rendering
	std::vector<float> vertices();
	std::vector<unsigned int> indices();
	int show(float dt, int n);

};

template int draw<Fluid_Euler>(Fluid_Euler* obj, float dt, int n);

#endif
