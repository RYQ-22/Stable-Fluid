#ifndef FLUID_EULER_H
#define FLUID_EULER_H

#include "draw.h"
#include "vector3.h"
#include "util.h"

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <functional>

#include <vector>
#include <unordered_map>
#include <queue>
#include <thread>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#define PI 3.1415926f
#define G_ 9.8f * 3
#define Gamma_ 0.1f

class Grid {
private:
	float l_;
	int N1_, N2_, N3_;
	int N_fluid_ = 0;
	// position of the grid's center
	Vec3<float> x_, y_, z_;
	// velocity at the grid's edge
	Vec3<float> Vx_, Vy_, Vz_;
	Vec3<float> Vx_temp, Vy_temp, Vz_temp;
	// level set
	Vec3<float> phi_;
	Vec3<int> particles_num_;
	// boundary
	//Vec3<int> solid_;
	Vec3<float> solid_phi_;
	Vec3<int> idx_;
	Vec3<int> Vx_valid_;
	Vec3<int> Vy_valid_;
	Vec3<int> Vz_valid_;
	// matrix
	Vec3<int> Adiag, Aplusi, Aplusj, Aplusk;
	Vec3<float> d, p;
	// z: auxiliary vetor
	// s: search vector
	Vec3<float> precon, q, z, s;

public:
	Grid();
	Grid(int N1, int N2, int N3, float l, std::vector<float> phi, std::vector<float> solid_phi);
	float Vx_ave(float i, float j, float k);
	float Vy_ave(float i, float j, float k);
	float Vz_ave(float i, float j, float k);
	float phi_ave(float i, float j, float k);
	float solid_phi_ave(float i, float j, float k);
	float solid_phi_ave(const glm::vec3& pos);


	friend class Fluid_Euler;
};

class Fluid_Euler {
private:
	float nu_ = 0.1f;
	float l_;	
	Grid grid_;

	float particle_radius_;
	std::vector<glm::vec3> particles_;
	std::vector<int> grid_idx_;
	std::vector<float> test_val_;

	// timing
	float deltaTime = 0.0f;
	float lastTime = 0.0f;
	float deltaTime_tot = 0.0f;
	// calculate frame rate
	float frameRate = 0;
	int count = 0;

public:
	Fluid_Euler(int N1, int N2, int N3, float l, std::vector<float> phi, std::vector<float> solid_phi);

	bool valid(int i, int j, int k);// test whether (i, j, k) is liquid

	void update(float dt);

	glm::vec3 trace_rk2(const glm::vec3& position, float dt);
	glm::vec3 get_velocity(const glm::vec3& position);

	// update fluid
	void add_force(float dt);// 1. add force
	void advect(float dt);// 2. advect
	void advect_particles(float dt);
	void project();// 3. project (matrix-free)
	void project0();
	void solve(int max_iterations);// solve: A * p = d;
	void applyPrecon();

	template <typename T>
	T dot(const Vec3<T>& v1, const Vec3<T>& v2);
	template <typename T>
	void applyA(const Vec3<T>& x, Vec3<T>& ans);

	// some details
	float compute_dt();
	void compute_phi(float dt);	
	void extrapolate();
	void constrain();

	// rendering
	std::vector<float> vertices();
	std::vector<unsigned int> indices();
	int show(float dt, int n);
	void run(float dt, int n);// only update, no ui
	void outputPLY(int fps, int t, float dt, int n);
	void simple_render(std::string output_file);
	void pbrt_render(std::string output_file);
};

template int draw<Fluid_Euler>(Fluid_Euler* obj, float dt, int n);

#endif
