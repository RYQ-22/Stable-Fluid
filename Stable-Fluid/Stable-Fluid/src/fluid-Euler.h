#ifndef FLUID_EULER_H
#define FLUID_EULER_H

#include "draw.h"
#include "field.h"
#include "util.h"
#include "solver.h"
#include "rigidbody.h"

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <math.h>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <functional>
#include <omp.h>

#include <vector>
#include <thread>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#ifndef PI
#define PI 3.1415926f
#endif
#ifndef G_
#define G_ 9.8f * 5
#endif
#ifndef Gamma_
#define Gamma_ 0.4f
#endif

class Grid {
private:
	float l_;
	int N1_, N2_, N3_;
	int N_fluid_ = 0;
	// position of the grid's center
	Field3<float> x_, y_, z_;
	// velocity at the grid's edge
	Field3<float> Vx_, Vy_, Vz_;
	Field3<float> new_Vx_, new_Vy_, new_Vz_;
	Field3<float> new_Vx_aux_, new_Vy_aux_, new_Vz_aux_;
	Field3<float> Vx_reflected_, Vy_reflected_, Vz_reflected_;
	// level set
	Field3<float> phi_;	
	Field3<float> new_phi_;
	// boundary
	//Vec3<int> solid_;
	Field3<float> solid_phi_;
	Field3<int> idx_; // only for project0
	Field3<int> Vx_valid_, Vy_valid_, Vz_valid_;
	Field3<int> old_Vx_valid_, old_Vy_valid_, old_Vz_valid_;

	// matrix
	Field3<float> Adiag, Aplusi, Aplusj, Aplusk;
	Field3<float> d, p;
	// z: auxiliary vetor
	// s: search vector
	Field3<float> precon, q, z, s;	

	// rigid body
	Field3<float> rigidBody_phi_;
	Field3<float> new_solid_phi_;

	Field3<float> J_x, J_y, J_z, Jrot_x, Jrot_y, Jrot_z;
	Field3<float> apply_J;

public:
	Grid();
	Grid(int N1, int N2, int N3, float l, std::vector<float> phi, std::vector<float> solid_phi);
	float Vx_ave(float i, float j, float k);
	float Vy_ave(float i, float j, float k);
	float Vz_ave(float i, float j, float k);
	float phi_ave(float i, float j, float k);
	float solid_phi_ave(float i, float j, float k);
	float solid_phi_ave(const glm::vec3& pos);
	float rigidBody_phi_ave(float i, float j, float k);
	float rigidBody_phi_ave(const glm::vec3& pos);

	friend class Fluid_Euler;
};

class Fluid_Euler {
private:
	float rho_ = 1.0f;// density of liquid
	float nu_ = 0.1f;
	float surfaceTensorCoefficient_ = 0.05f;
	float l_;	
	Grid grid_;

	// setting
	bool use_mc_ = false;
	bool use_reflection_ = false;
	bool add_surfaceTension_ = false;
	bool use_fixedOperation_ = false;

	// particles
	float particle_radius_;
	std::vector<glm::vec3> particles_;

	// timing
	double deltaTime_ = 0.0;
	double lastTime_ = 0.0;
	double deltaTime_tot_ = 0.0;
	// calculate frame rate
	double frameRate_ = 0;
	int count_ = 0;

	// solid-fluid coupling
	bool add_rigidBody_ = false;
	RigidBody rigidBody_;
	std::vector<glm::vec3> rigidBody_particles_;
	
public:
	Fluid_Euler(int N1, int N2, int N3, float l, std::vector<float> phi, std::vector<float> solid_phi);
	void set_velocity(glm::vec3 v);

	bool valid(int i, int j, int k);// test whether (i, j, k) is liquid

	template <typename T>
	T dot(const Field3<T>& v1, const Field3<T>& v2);	
	void applyA(const Field3<float>& x, Field3<float>& ans);	
	Field3<float> applyA(const Field3<float>& x);
	
	// update fluid
	void update(float dt);
	// 1. add force
	void add_force(float dt);
	// 2. advect
	// id = 0: v_
	// id = 1: v_reflected
	void advect(float dt, int field_id = 0);
	void advect_particles(float dt);

	glm::vec3 get_velocity(const glm::vec3& position);
	glm::vec3 trace_rk2(const glm::vec3& position, float dt);
	void semi_Lagrangian(const Field3<float>& field, Field3<float>& field_new, float dt);
	void mac_Cormack(const Field3<float>& field, Field3<float>& new_field, Field3<float>& new_field_aux, float dt);
	// 3. project (matrix-free)
	void project();
	void apply_pressure_liquid(Field3<float>& vx, Field3<float>& vy, Field3<float>& vz);
	void project_new();
	void project0();// use eigen
	void solve(int max_iterations);// solve: A * p = d;
	void applyPrecon();

	// some details
	float compute_dt();
	void compute_phi(float dt);	
	void extrapolate();
	void extrapolate(Field3<float>& v, Field3<float> v_new, Field3<int>& valid, Field3<int>& valid_old);
	void constrain();

	// solid-fluid coupling
	void add_rigidBody(// add rigidBody from .obj
		std::string obj_path,
		float scale, 
		glm::vec3 translate,
		float M,
		glm::mat3 I_ref,
		glm::vec3 c0,
		glm::quat rotationQuaternion = glm::quat(1.0f, 0.0f, 0.0f, 0.0f),
		glm::vec3 Vc = glm::vec3(0.0f, 0.0f, 0.0f),
		glm::vec3 omega = glm::vec3(0.0f, 0.0f, 0.0f));
	void update_rigidBody(float dt);
	void set_fixedOperation();
	float computeVolume(int i, int j, int k, const Field3<float>& phi, int id);
	void apply_pressure_rigidBody();

	// rendering
	std::vector<float> vertices();
	std::vector<unsigned int> indices();
	int show(float dt, int n);
	void run(float dt, int n);// only update, no ui
	void outputPLY(int fps, int t, float dt, int n);
	void outputXYZ(int fps, int t, float dt, int n);	
	void simple_render(std::string output_file, int fps = 30);
	void pbrt_render(std::string output_file, int fps, int t);
};


template int draw<Fluid_Euler>(Fluid_Euler* obj, float dt, int n);

class CleanupHelper {
public:
	CleanupHelper(const std::string& folderPath) : folderPath_(folderPath) {}

	~CleanupHelper() {
		try {
			for (const auto& entry : std::filesystem::directory_iterator(folderPath_)) {
				std::filesystem::remove(entry.path());
			}
		}
		catch (const std::filesystem::filesystem_error& err) {
			std::cerr << "Error during cleanup: " << err.what() << std::endl;
		}
	}

private:
	std::string folderPath_;
};

#endif
