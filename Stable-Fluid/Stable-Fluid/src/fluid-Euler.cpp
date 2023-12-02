#include"fluid-Euler.h"

Grid::Grid() {}

Grid::Grid(int N1, int N2, int N3, float l, std::vector<float> phi, std::vector<float> solid_phi) : N1_(N1), N2_(N2), N3_(N3), l_(l) {
	if (phi.size() != N1 * N2 * N3) {
		std::cout << "init level set error!" << std::endl;
		abort();
	}
	if (solid_phi.size() != N1 * N2 * N3) {
		std::cout << "init solid set error!" << std::endl;
		abort();
	}

	x_.resize(N1, N2, N3);
	y_.resize(N1, N2, N3);
	z_.resize(N1, N2, N3);

	Vx_.resize((N1 + 1), N2, N3);
	Vy_.resize(N1, (N2 + 1), N3);
	Vz_.resize(N1, N2, (N3 + 1));
	Vx_temp.resize((N1 + 1), N2, N3);
	Vy_temp.resize(N1, (N2 + 1), N3);
	Vz_temp.resize(N1, N2, (N3 + 1));

	idx_.resize(N1, N2, N3);

	Vx_valid_.resize(N1 + 1, N2, N3);
	Vy_valid_.resize(N1, N2 + 1, N3);
	Vz_valid_.resize(N1, N2, N3 + 1);

	Adiag.resize(N1, N2, N3);
	Aplusi.resize(N1, N2, N3);
	Aplusj.resize(N1, N2, N3);
	Aplusk.resize(N1, N2, N3);
	d.resize(N1, N2, N3);
	p.resize(N1, N2, N3);
	precon.resize(N1, N2, N3);
	q.resize(N1, N2, N3);
	z.resize(N1, N2, N3);
	s.resize(N1, N2, N3);

	phi_ = Vec3<float>(N1, N2, N3, phi);
	solid_phi_ = Vec3<float>(N1, N2, N3, solid_phi);

	for (int i = 0; i < N1; i++) for (int j = 0; j < N2; j++) for (int k = 0; k < N3; k++) {
		// init position
		x_(i, j, k) = i * l_;
		y_(i, j, k) = j * l_;
		z_(i, j, k) = k * l_;
		// init N_fluid_
		if (solid_phi_(i, j, k) > 0 && phi_(i, j, k) < 0) {
			N_fluid_++;
		}
	}		
}

float Grid::Vx_ave(float i, float j, float k) {
	return interpolate_value(i, j, k, Vx_);
}

float Grid::Vy_ave(float i, float j, float k) {
	return interpolate_value(i, j, k, Vy_);
}

float Grid::Vz_ave(float i, float j, float k) {
	return interpolate_value(i, j, k, Vz_);
}

float Grid::phi_ave(float i, float j, float k) {
	return interpolate_value(i, j, k, phi_);
}

float Grid::solid_phi_ave(float i, float j, float k) {
	return interpolate_value(i, j, k, solid_phi_);
}

float Grid::solid_phi_ave(glm::vec3 pos) {
	return interpolate_value(pos.x, pos.y, pos.z, solid_phi_);
}

Fluid_Euler::Fluid_Euler(int N1, int N2, int N3, float l, std::vector<float> phi, std::vector<float> solid_phi) : l_(l) {
	grid_ = Grid(N1, N2, N3, l, phi, solid_phi);

	// make the particles large enough so they always appear on the grid
	particle_radius_ = (float)(l_ * 1.01 * sqrt(3.0) / 2.0);
	// init particles
	int seed = 0;
	for (int n = 0; n < 4; n++) {
		for (int i = 0; i < N1; i++) for (int j = 0; j < N2; j++) for (int k = 0; k < N3; k++) {
			float a = randhashf(seed++); float b = randhashf(seed++); float c = randhashf(seed++);
			float x = (float)i + a, y = (float)j + b, z = (float)k + c;
			if (grid_.phi_ave(x, y, z) <= -particle_radius_) {				
				if (grid_.solid_phi_ave(x, y, z) > 0)
					particles_.push_back(glm::vec3(x * l, y * l, z * l));
			}
		}
	}	
}

// Apply RK2 to advect a point in the domain.
glm::vec3 Fluid_Euler::trace_rk2(const glm::vec3& position, float dt) {
	glm::vec3 input = position;
	glm::vec3 velocity = get_velocity(input);
	velocity = get_velocity(input + 0.5f * dt * velocity);
	input += dt * velocity;
	return input;
}

glm::vec3 Fluid_Euler::get_velocity(const glm::vec3& position) {
	float vx_ave = grid_.Vx_ave(position.x / l_ + 0.5f, position.y / l_, position.z / l_);
	float vy_ave = grid_.Vy_ave(position.x / l_, position.y / l_ + 0.5f, position.z / l_);
	float vz_ave = grid_.Vz_ave(position.x / l_, position.y / l_, position.z / l_ + 0.5f);
	return glm::vec3(vx_ave, vy_ave, vz_ave);
}

bool Fluid_Euler::valid(int i, int j, int k) {// fluid, and not in solid
	if (i >= 0 && i < grid_.N1_ && j >= 0 && j < grid_.N2_ && k >= 0 && k < grid_.N3_ && grid_.solid_phi_(i, j, k) > 0 && grid_.phi_(i, j, k) < 0.0f) {
		return true;
	}
	return false;
}

template <typename T>
T Fluid_Euler::dot(const Vec3<T>& v1, const Vec3<T>& v2) {
	if (v1.N1() != v2.N1() || v1.N2() != v2.N2() || v1.N3() != v2.N3()) {
		std::cerr << "2 Vec3 have diffrent size!" << std::endl;
	}
	T ans = 0;
	int n1 = v1.N1();
	int n2 = v1.N2();
	int n3 = v1.N3();
	for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++) {
		if (valid(i, j, k))
			ans += v1(i, j, k) * v2(i, j, k);
	}
	return ans;
}

template <typename T>
void Fluid_Euler::applyA(const Vec3<T>& x, Vec3<T>& ans) {
	for (int i = 0; i < x.N1(); i++) for (int j = 0; j < x.N2(); j++) for (int k = 0; k < x.N3(); k++) {
		ans(i, j, k) = grid_.Adiag(i, j, k) * x(i, j, k);
		if (valid(i + 1, j, k))
			ans(i, j, k) += grid_.Aplusi(i, j, k) * x(i + 1, j, k);
		if (valid(i, j + 1, k))
			ans(i, j, k) += grid_.Aplusj(i, j, k) * x(i, j + 1, k);
		if (valid(i, j, k + 1))
			ans(i, j, k) += grid_.Aplusk(i, j, k) * x(i, j, k + 1);
		if (valid(i - 1, j, k))
			ans(i, j, k) += grid_.Aplusi(i - 1, j, k) * x(i - 1, j, k);
		if (valid(i, j - 1, k))
			ans(i, j, k) += grid_.Aplusj(i, j - 1, k) * x(i, j - 1, k);
		if (valid(i, j, k - 1))
			ans(i, j, k) += grid_.Aplusk(i, j, k - 1) * x(i, j, k - 1);
	}
	return;
}

void Fluid_Euler::add_force(float dt) {
	for (int i = 0; i < grid_.N1_ + 1; i++) for (int j = 0; j < grid_.N2_ + 1; j++) for (int k = 0; k < grid_.N3_ + 1; k++) {
		// Vy
		if (i < grid_.N1_ && k < grid_.N3_) {
			if ((j < grid_.N2_ && valid(i, j, k)) || (j != 0 && valid(i, j - 1, k))) 
				grid_.Vy_(i, j, k) += -G_ * dt;			
		}
	}
	return;
}


void Fluid_Euler::advect_particles(float dt) {
	for (int i = 0; i < particles_.size(); i++) {
		particles_[i] = trace_rk2(particles_[i], dt);
		// boundary
		float phi_val = grid_.solid_phi_ave(particles_[i]);
		if (phi_val < 0) {
			glm::vec3 grad;
			interpolate_gradient(grad, particles_[i] / l_, grid_.solid_phi_);
			grad = glm::normalize(grad);			
			//particles_[i] -= phi_val * grad;
		}
	}
}


void Fluid_Euler::advect(float dt) {		
	for (int i = 0; i < grid_.N1_ + 1; i++) for (int j = 0; j < grid_.N2_ + 1; j++) for (int k = 0; k < grid_.N3_ + 1; k++) {
		// Vx
		if (j < grid_.N2_ && k < grid_.N3_) {
			glm::vec3 pos((i - 0.5f) * l_, j * l_, k * l_);// from vx coordinate to center coordinate
			pos = trace_rk2(pos, -dt);
			grid_.Vx_temp(i, j, k) = get_velocity(pos).x;
		}			
		// Vy
		if (i < grid_.N1_ && k < grid_.N3_) {
			glm::vec3 pos(i * l_, (j - 0.5f) * l_, k * l_);
			pos = trace_rk2(pos, -dt);
			grid_.Vy_temp(i, j, k) = get_velocity(pos).y;
		}			
		// Vz
		if (i < grid_.N1_ && j < grid_.N2_) {
			glm::vec3 pos(i * l_, j * l_, (k - 0.5f) * l_);
			pos = trace_rk2(pos, -dt);
			grid_.Vz_temp(i, j, k) = get_velocity(pos).z;
		}			
	}
			
	grid_.Vx_ = grid_.Vx_temp;
	grid_.Vy_ = grid_.Vy_temp;
	grid_.Vz_ = grid_.Vz_temp;
	return;
}

// use eigen
void Fluid_Euler::project0() {
	// compute number of fluid grids
	grid_.N_fluid_ = 0;
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
		if (valid(i, j, k)) {
			grid_.idx_(i, j, k) = grid_.N_fluid_;
			grid_.N_fluid_++;
		}
	}

	Eigen::VectorXf p(grid_.N_fluid_);
	Eigen::VectorXf d(grid_.N_fluid_);
	Eigen::SparseMatrix<float> A(grid_.N_fluid_, grid_.N_fluid_);

	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
		if (valid(i, j, k)) {
			int neighbor = 0;
			float di = 0.0f;
			int index = grid_.idx_(i, j, k);
			for (int i1 = -1; i1 <= 1; i1 += 2) {
				if (grid_.solid_phi_(i + i1, j, k) > 0) {
					neighbor++;
				}
				if (grid_.solid_phi_(i, j + i1, k) > 0) {
					neighbor++;
				}
				if (grid_.solid_phi_(i, j, k + i1) > 0) {
					neighbor++;
				}				
				if (valid(i + i1, j, k)) {// fluid neighbor
					A.insert(index, grid_.idx_(i + i1, j, k)) = -1.0f;
				}
				if (valid(i, j + i1, k)) {
					A.insert(index, grid_.idx_(i, j + i1, k)) = -1.0f;
				}
				if (valid(i, j, k + i1)) {
					A.insert(index, grid_.idx_(i, j, k + i1)) = -1.0f;
				}
				// calculate divergence of velocity
				if (i1 == 1) {
					di += grid_.solid_phi_(i + 1, j, k) > 0 ? -grid_.Vx_(i + 1, j, k) : 0;
					di += grid_.solid_phi_(i, j + 1, k) > 0 ? -grid_.Vy_(i, j + 1, k) : 0;
					di += grid_.solid_phi_(i, j, k + 1) > 0 ? -grid_.Vz_(i, j, k + 1) : 0;
				}
				if (i1 == -1) {
					di += grid_.solid_phi_(i - 1, j, k) > 0 ? grid_.Vx_(i, j, k) : 0;
					di += grid_.solid_phi_(i, j - 1, k) > 0 ? grid_.Vy_(i, j, k) : 0;
					di += grid_.solid_phi_(i, j, k - 1) > 0 ? grid_.Vz_(i, j, k) : 0;
				}
			}
			A.insert(index, index) = neighbor;
			d(index) = di;
		}
	}

	A.makeCompressed();

	// CG
	Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> cg;
	//cg.setMaxIterations(10);
	cg.compute(A);
	p = cg.solve(d);
	// update Vx, Vy, Vz
	//std::cout << grid_.Vy_(10,10,10) << std::endl;

	for (int i = 0; i < grid_.N1_ + 1; i++) for (int j = 0; j < grid_.N2_ + 1; j++) for (int k = 0; k < grid_.N3_ + 1; k++) {
		// Vx
		if (j < grid_.N2_ && k < grid_.N3_) {
			if ((i < grid_.N1_ && valid(i, j, k)) || (i != 0 && valid(i - 1, j, k))) {
				if (!(i < grid_.N1_ && grid_.solid_phi_(i, j, k) < 0) && !(i != 0 && grid_.solid_phi_(i - 1, j, k) < 0)) {// not solid
					if (i < grid_.N1_ && grid_.phi_(i, j, k) < 0.0f)// p+
						grid_.Vx_(i, j, k) += -p(grid_.idx_(i, j, k));
					if (i != 0 && grid_.phi_(i - 1, j, k) < 0.0f)// p-
						grid_.Vx_(i, j, k) += p(grid_.idx_(i - 1, j, k));
				}
			}
		}
		// Vy
		if (i < grid_.N1_ && k < grid_.N3_) {
			if ((j < grid_.N2_ && valid(i, j, k)) || (j != 0 && valid(i, j - 1, k))) {
				if (!(j < grid_.N2_ && grid_.solid_phi_(i, j, k) < 0) && !(j != 0 && grid_.solid_phi_(i, j - 1, k) < 0)) {// not solid
					if (j < grid_.N2_ && grid_.phi_(i, j, k) < 0.0f)// p+
						grid_.Vy_(i, j, k) += -p(grid_.idx_(i, j, k));
					if (j != 0 && grid_.phi_(i, j - 1, k) < 0.0f)// p-
						grid_.Vy_(i, j, k) += p(grid_.idx_(i, j - 1, k));
				}
			}
		}
		// Vz
		if (i < grid_.N1_ && j < grid_.N2_) {
			if ((k < grid_.N3_ && valid(i, j, k)) || (k != 0 && valid(i, j, k - 1))) {
				if (!(k < grid_.N3_ && grid_.solid_phi_(i, j, k) < 0) && !(k != 0 && grid_.solid_phi_(i, j, k - 1) < 0)) {// not solid
					if (k < grid_.N3_ && grid_.phi_(i, j, k) < 0.0f)// p+
						grid_.Vz_(i, j, k) += -p(grid_.idx_(i, j, k));
					if (k != 0 && grid_.phi_(i, j, k - 1) < 0.0f)// p-
						grid_.Vz_(i, j, k) += p(grid_.idx_(i, j, k - 1));
				}
			}
		}
	}	
}

// matrix-free CG
void Fluid_Euler::project() {
	// compute A, d
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {		
		if (valid(i, j, k)) {// for liquid
			grid_.Adiag(i, j, k) = 0;
			grid_.Aplusi(i, j, k) = 0;
			grid_.Aplusj(i, j, k) = 0;
			grid_.Aplusk(i, j, k) = 0;
			grid_.d(i, j, k) = 0;
			for (int i1 = -1; i1 <= 1; i1 += 2) {
				if (grid_.solid_phi_(i + i1, j, k) > 0) {
					grid_.Adiag(i, j, k)++;
				}
				if (grid_.solid_phi_(i, j + i1, k) > 0) {
					grid_.Adiag(i, j, k)++;
				}
				if (grid_.solid_phi_(i, j, k + i1) > 0) {
					grid_.Adiag(i, j, k)++;
				}				
				if(i1 == 1) {
					grid_.d(i, j, k) += grid_.solid_phi_(i + 1, j, k) > 0 ? -grid_.Vx_(i + 1, j, k) : 0;
					grid_.d(i, j, k) += grid_.solid_phi_(i, j + 1, k) > 0 ? -grid_.Vy_(i, j + 1, k) : 0;
					grid_.d(i, j, k) += grid_.solid_phi_(i, j, k + 1) > 0 ? -grid_.Vz_(i, j, k + 1) : 0;
				}
				if (i1 == -1) {
					grid_.d(i, j, k) += grid_.solid_phi_(i - 1, j, k) > 0 ? grid_.Vx_(i, j, k) : 0;
					grid_.d(i, j, k) += grid_.solid_phi_(i, j - 1, k) > 0 ? grid_.Vy_(i, j, k) : 0;
					grid_.d(i, j, k) += grid_.solid_phi_(i, j, k - 1) > 0 ? grid_.Vz_(i, j, k) : 0;
				}
			}
			if (valid(i + 1, j, k)) {
				grid_.Aplusi(i, j, k) = -1.0f;
			}
			if (valid(i, j + 1, k)) {
				grid_.Aplusj(i, j, k) = -1.0f;
			}
			if (valid(i, j, k + 1)) {
				grid_.Aplusk(i, j, k) = -1.0f;
			}
		}
	}
	// solve: A * p = d
	solve(20);
	// update Vx, Vy, Vz
	for (int i = 0; i < grid_.N1_ + 1; i++) for (int j = 0; j < grid_.N2_ + 1; j++) for (int k = 0; k < grid_.N3_ + 1; k++) {
		// Vx
		if (j < grid_.N2_ && k < grid_.N3_) {
			if ((i < grid_.N1_ && valid(i, j, k)) || (i != 0 && valid(i - 1, j, k))) {// liquid neighbor
				if (!(i < grid_.N1_ && grid_.solid_phi_(i, j, k) < 0) && !(i != 0 && grid_.solid_phi_(i - 1, j, k) < 0)) {// not solid
					if (i < grid_.N1_ && grid_.phi_(i, j, k) < 0.0f)// p+
						grid_.Vx_(i, j, k) += -grid_.p(i, j, k);
					if (i != 0 && grid_.phi_(i - 1, j, k) < 0.0f)// p-
						grid_.Vx_(i, j, k) += grid_.p(i - 1, j, k);
				}
			}
		}
		// Vy
		if (i < grid_.N1_ && k < grid_.N3_) {
			if ((j < grid_.N2_ && valid(i, j, k)) || (j != 0 && valid(i, j - 1, k))) {
				if (!(j < grid_.N2_ && grid_.solid_phi_(i, j, k) < 0) && !(j != 0 && grid_.solid_phi_(i, j - 1, k) < 0)) {// not solid
					if (j < grid_.N2_ && grid_.phi_(i, j, k) < 0.0f)// p+
						grid_.Vy_(i, j, k) += -grid_.p(i, j, k);
					if (j != 0 && grid_.phi_(i, j - 1, k) < 0.0f)// p-
						grid_.Vy_(i, j, k) += grid_.p(i, j - 1, k);
				}
			}
		}
		// Vz
		if (i < grid_.N1_ && j < grid_.N2_) {
			if ((k < grid_.N3_ && valid(i, j, k)) || (k != 0 && valid(i, j, k - 1))) {
				if (!(k < grid_.N3_ && grid_.solid_phi_(i, j, k) < 0) && !(k != 0 && grid_.solid_phi_(i, j, k - 1) < 0)) {// not solid
					if (k < grid_.N3_ && grid_.phi_(i, j, k) < 0.0f)// p+
						grid_.Vz_(i, j, k) += -grid_.p(i, j, k);
					if (k != 0 && grid_.phi_(i, j, k - 1) < 0.0f)// p-
						grid_.Vz_(i, j, k) += grid_.p(i, j, k - 1);
				}
			}
		}
	}
}

void Fluid_Euler::applyPrecon() {
	// compute precon
	float tau = 0.97f;
	float rho = 0.25f;
	float e = 0.0f;
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
		if (valid(i, j, k)) {
			e = grid_.Adiag(i, j, k);
			e += -pow(grid_.Aplusi(i - 1, j, k) * grid_.precon(i - 1, j, k), 2) - tau * grid_.Aplusi(i - 1, j, k) * (grid_.Aplusj(i - 1, j, k) + grid_.Aplusk(i - 1, j, k)) * pow(grid_.precon(i - 1, j, k), 2);
			e += -pow(grid_.Aplusj(i, j - 1, k) * grid_.precon(i, j - 1, k), 2) - tau * grid_.Aplusj(i, j - 1, k) * (grid_.Aplusi(i, j - 1, k) + grid_.Aplusk(i, j - 1, k)) * pow(grid_.precon(i, j - 1, k), 2);
			e += -pow(grid_.Aplusk(i, j, k - 1) * grid_.precon(i, j, k - 1), 2) - tau * grid_.Aplusk(i, j, k - 1) * (grid_.Aplusi(i, j, k - 1) + grid_.Aplusj(i, j, k - 1)) * pow(grid_.precon(i, j, k - 1), 2);
			if (e < rho * grid_.Adiag(i, j, k))
				e = grid_.Adiag(i, j, k);
			grid_.precon(i, j, k) = 1 / sqrt(e + 1e-10f);
		}
		else {
			grid_.precon(i, j, k) = 0;
		}
	}
	// solve Lq = d
	float t = 0.0f;
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
		if (valid(i, j, k)) {
			t = grid_.d(i, j, k);
			if (valid(i - 1, j, k))
				t += -grid_.Aplusi(i - 1, j, k) * grid_.precon(i - 1, j, k) * grid_.q(i - 1, j, k);
			if (valid(i, j - 1, k))
				t += -grid_.Aplusj(i, j - 1, k) * grid_.precon(i, j - 1, k) * grid_.q(i, j - 1, k);
			if (valid(i, j, k - 1))
				t += -grid_.Aplusk(i, j, k - 1) * grid_.precon(i, j, k - 1) * grid_.q(i, j, k - 1);
			grid_.q(i, j, k) = t * grid_.precon(i, j, k);
		}
	}
	//solve L^T z = q
	for (int i = grid_.N1_ - 1; i >= 0; i--) for (int j = grid_.N2_ - 1; j >= 0; j--) for (int k = grid_.N3_ - 1; k >= 0; k--) {
		if (valid(i, j, k)) {
			t = grid_.q(i, j, k);
			if (valid(i + 1, j, k))
				t += -grid_.Aplusi(i, j, k) * grid_.precon(i, j, k) * grid_.z(i + 1, j, k);
			if (valid(i, j + 1, k))
				t += -grid_.Aplusj(i, j, k) * grid_.precon(i, j, k) * grid_.z(i, j + 1, k);
			if (valid(i, j, k + 1))
				t += -grid_.Aplusk(i, j, k) * grid_.precon(i, j, k) * grid_.z(i, j, k + 1);
			grid_.z(i, j, k) = t * grid_.precon(i, j, k);
		}
	}
}

void Fluid_Euler::solve(int max_iterations) {
	// set p = 0
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {		
		grid_.p(i, j, k) = 0.0f;	
	}
	if (grid_.d.max() == 0)
		return;
	// PCG
	applyPrecon();
	grid_.s = grid_.z;
	float sigma = dot(grid_.z, grid_.d);
	float sigma_new;
	float alpha;
	for (int i = 0; i < max_iterations; i++) {
		applyA(grid_.s, grid_.z);
		alpha = sigma / dot(grid_.z, grid_.s);
		grid_.p += grid_.s * alpha;
		grid_.d -= grid_.z * alpha;
		if (grid_.d.max() <= 1e-6f)
			return;
		applyPrecon();		
		sigma_new = dot(grid_.z, grid_.d);
		alpha = sigma_new / sigma;
		sigma = sigma_new;
		grid_.s = grid_.z + grid_.s * alpha;
	}
	return;
}

void Fluid_Euler::compute_phi(float dt) {
	grid_.phi_.assign(3 * l_);
	for (int p = 0; p < particles_.size(); p++) {
		glm::ivec3 cell_ind(particles_[p] / l_);
		for (int k = max(0, cell_ind[2] - 1); k <= min(cell_ind[2] + 1, grid_.N3_ - 1); ++k) {
			for (int j = max(0, cell_ind[1] - 1); j <= min(cell_ind[1] + 1, grid_.N2_ - 1); ++j) {
				for (int i = max(0, cell_ind[0] - 1); i <= min(cell_ind[0] + 1, grid_.N1_ - 1); ++i) {
					glm::vec3 sample_pos((i ) * l_, (j ) * l_, (k ) * l_);					
					float test_val = glm::length(sample_pos-particles_[p]) - particle_radius_;
					if (test_val < grid_.phi_(i, j, k))
						grid_.phi_(i, j, k) = test_val;
				}
			}
		}
	}
	
	// extend phi into solid
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
		if (grid_.phi_(i, j, k) < 0.5f * l_ && grid_.solid_phi_(i,j,k) < 0) {
			grid_.phi_(i, j, k) = -0.5f * l_;
		}
	}
	
}

void Fluid_Euler::extrapolate() {
	// reset valid_
	for (int i = 0; i < grid_.N1_ + 1; i++) for (int j = 0; j < grid_.N2_ + 1; j++) for (int k = 0; k < grid_.N3_ + 1; k++) {
		if (j < grid_.N2_ && k < grid_.N3_) {
			if ((i < grid_.N1_ && grid_.solid_phi_(i, j, k) > 0 && grid_.phi_(i, j, k) < 0.0f) || (i - 1 >= 0 && grid_.solid_phi_(i - 1, j, k) > 0 && grid_.phi_(i - 1, j, k) < 0.0f)) {
				grid_.Vx_valid_(i, j, k) = 1;
			}
			else {
				grid_.Vx_valid_(i, j, k) = 0;
			}
		}
		if (i < grid_.N1_ && k < grid_.N3_) {
			if ((j < grid_.N2_ && grid_.solid_phi_(i, j, k) > 0 && grid_.phi_(i, j, k) < 0.0f) || (j - 1 >= 0 && grid_.solid_phi_(i, j - 1, k) > 0 && grid_.phi_(i, j - 1, k) < 0.0f)) {
				grid_.Vy_valid_(i, j, k) = 1;
			}
			else {
				grid_.Vy_valid_(i, j, k) = 0;
			}
		}
		if (i < grid_.N1_ && j < grid_.N2_) {
			if ((k < grid_.N3_ && grid_.solid_phi_(i, j, k) > 0 && grid_.phi_(i, j, k) < 0.0f) || (k - 1 >= 0 && grid_.solid_phi_(i, j, k - 1) < 0 && grid_.phi_(i, j, k - 1) < 0.0f)) {
				grid_.Vz_valid_(i, j, k) = 1;
			}
			else {
				grid_.Vz_valid_(i, j, k) = 0;
			}
		}
	}
	// Apply several iterations of a very simple propagation of valid velocity data in all directions
	for (int num = 0; num < 15; num++) {
		for (int i = 0; i < grid_.N1_ + 1; i++) for (int j = 0; j < grid_.N2_ + 1; j++) for (int k = 0; k < grid_.N3_ + 1; k++) {
			float vx = 0, vy = 0, vz = 0;
			int countx = 0, county = 0, countz = 0;
			// Vx
			if (j < grid_.N2_ && k < grid_.N3_) {
				if (!grid_.Vx_valid_(i, j, k)) {
					if (grid_.Vx_valid_(i - 1, j, k)) {
						vx += grid_.Vx_(i - 1, j, k);
						countx++;
					}
					if (grid_.Vx_valid_(i + 1, j, k)) {
						vx += grid_.Vx_(i + 1, j, k);
						countx++;
					}
					if (grid_.Vx_valid_(i, j - 1, k)) {
						vx += grid_.Vx_(i, j - 1, k);
						countx++;
					}
					if (grid_.Vx_valid_(i, j + 1, k)) {
						vx += grid_.Vx_(i, j + 1, k);
						countx++;
					}
					if (grid_.Vx_valid_(i, j, k - 1)) {
						vx += grid_.Vx_(i, j, k - 1);
						countx++;
					}
					if (grid_.Vx_valid_(i, j, k + 1)) {
						vx += grid_.Vx_(i, j, k + 1);
						countx++;
					}
					if (countx > 0) {
						grid_.Vx_(i, j, k) = vx / countx;
						grid_.Vx_valid_(i, j, k) = 1;
					}
				}
			}
			// Vy
			if (i < grid_.N1_ && k < grid_.N3_) {
				if (!grid_.Vy_valid_(i, j, k)) {
					if (grid_.Vy_valid_(i - 1, j, k)) {
						vy += grid_.Vy_(i - 1, j, k);
						county++;
					}
					if (grid_.Vy_valid_(i + 1, j, k)) {
						vy += grid_.Vy_(i + 1, j, k);
						county++;
					}
					if (grid_.Vy_valid_(i, j - 1, k)) {
						vy += grid_.Vy_(i, j - 1, k);
						county++;
					}
					if (grid_.Vy_valid_(i, j + 1, k)) {
						vy += grid_.Vy_(i, j + 1, k);
						county++;
					}
					if (grid_.Vy_valid_(i, j, k - 1)) {
						vy += grid_.Vy_(i, j, k - 1);
						county++;
					}
					if (grid_.Vy_valid_(i, j, k + 1)) {
						vy += grid_.Vy_(i, j, k + 1);
						county++;
					}
					if (county > 0) {
						grid_.Vy_(i, j, k) = vy / county;
						grid_.Vy_valid_(i, j, k) = 1;
					}
				}
			}
			// Vz
			if (i < grid_.N1_ && j < grid_.N2_) {
				if (!grid_.Vz_valid_(i, j, k)) {
					if (grid_.Vz_valid_(i - 1, j, k)) {
						vz += grid_.Vz_(i - 1, j, k);
						countz++;
					}
					if (grid_.Vz_valid_(i + 1, j, k)) {
						vz += grid_.Vz_(i + 1, j, k);
						countz++;
					}
					if (grid_.Vz_valid_(i, j - 1, k)) {
						vz += grid_.Vz_(i, j - 1, k);
						countz++;
					}
					if (grid_.Vz_valid_(i, j + 1, k)) {
						vz += grid_.Vz_(i, j + 1, k);
						countz++;
					}
					if (grid_.Vz_valid_(i, j, k - 1)) {
						vz += grid_.Vz_(i, j, k - 1);
						countz++;
					}
					if (grid_.Vz_valid_(i, j, k + 1)) {
						vz += grid_.Vz_(i, j, k + 1);
						countz++;
					}
					if (countz > 0) {
						grid_.Vz_(i, j, k) = vz / countz;
						grid_.Vz_valid_(i, j, k) = 1;
					}
				}
			}
		}
	}
}


void Fluid_Euler::constrain() {
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
		if (valid(i, j, k)) {// for liquid
			if ((i + 1 < grid_.N1_ && grid_.solid_phi_(i + 1, j, k) < 0 && grid_.Vx_(i + 1, j, k) > 0))
				grid_.Vx_(i + 1, j, k) = 0;
			if (i > 0 && grid_.solid_phi_(i - 1, j, k) < 0 && grid_.Vx_(i, j, k) < 0)
				grid_.Vx_(i, j, k) = 0;
			if ((j + 1 < grid_.N2_ && grid_.solid_phi_(i, j + 1, k) < 0 && grid_.Vy_(i, j + 1, k) > 0))
				grid_.Vy_(i, j + 1, k) = 0;
			if (j > 0 && grid_.solid_phi_(i, j - 1, k) < 0 && grid_.Vy_(i, j, k) < 0)
				grid_.Vy_(i, j, k) = 0;
			if ((k + 1 < grid_.N3_ && grid_.solid_phi_(i, j, k) < 0 && grid_.Vz_(i, j, k + 1) > 0))
				grid_.Vz_(i, j, k + 1) = 0;
			if (k > 0 && grid_.solid_phi_(i, j, k - 1) < 0 && grid_.Vz_(i, j, k) < 0)
				grid_.Vz_(i, j, k) = 0;
		}
	}
}

float Fluid_Euler::compute_dt() {	
	float Vx_max = grid_.Vx_.max();
	float Vy_max = grid_.Vy_.max();
	float Vz_max = grid_.Vz_.max();

	return l_ / (std::max(Vx_max, std::max(Vy_max, Vz_max)) + 1e-6f);
}


void Fluid_Euler::update(float dt) {
	std::cout << grid_.N_fluid_ << std::endl;
	// 1. extrapolation
	extrapolate();
	// 2. update phi
	advect_particles(dt);
	compute_phi(dt);
	// 3. update v
	// advect
	advect(dt);
	// add force
	add_force(dt);
	// project A * p = d
	project();
	// constrain for fluid
	constrain();// cause unsymmetry
}


std::vector<float> Fluid_Euler::vertices() {
	std::vector<float> vertices;

	grid_.N_fluid_ = 0;
	for (int i0 = 0; i0 < grid_.N1_; i0++) for (int j0 = 0; j0 < grid_.N2_; j0++) for (int k0 = 0; k0 < grid_.N3_; k0++) {
		if (valid(i0, j0, k0)) {
			grid_.N_fluid_++;
			for (int i = -1; i <= 1; i += 2) for (int j = -1; j <= 1; j += 2) for (int k = -1; k <= 1; k += 2) {
				vertices.push_back(10 * (grid_.x_(i0, j0, k0) + 0.5f * l_ * i));
				vertices.push_back(10 * (grid_.y_(i0, j0, k0) + 0.5f * l_ * j));
				vertices.push_back(10 * (grid_.z_(i0, j0, k0) + 0.5f * l_ * k));
				vertices.push_back(0.67843f);
				vertices.push_back(0.84706f);
				vertices.push_back(0.90196f);
			}					
		}
	}

	/*
	//a plane y=0.0f
	vertices.push_back(-20.0f);
	vertices.push_back(0.0f);
	vertices.push_back(-20.0f);
	vertices.push_back(0.75f);
	vertices.push_back(0.75f);
	vertices.push_back(0.75f);

	vertices.push_back(-20.0f);
	vertices.push_back(0.0f);
	vertices.push_back(20.0f);
	vertices.push_back(0.25f);
	vertices.push_back(0.25f);
	vertices.push_back(0.25f);

	vertices.push_back(20.0f);
	vertices.push_back(0.0f);
	vertices.push_back(-20.0f);
	vertices.push_back(0.25f);
	vertices.push_back(0.25f);
	vertices.push_back(0.25f);

	vertices.push_back(20.0f);
	vertices.push_back(0.0f);
	vertices.push_back(20.0f);
	vertices.push_back(0.25f);
	vertices.push_back(0.25f);
	vertices.push_back(0.25f);
	*/
	return vertices;
}

std::vector<unsigned int> Fluid_Euler::indices() {
	std::vector<unsigned int> indices;

	for (int i = 0; i < grid_.N_fluid_; i++) {// may be bug
		indices.push_back(8 * i);
		indices.push_back(8 * i + 1);
		indices.push_back(8 * i + 5);

		indices.push_back(8 * i);
		indices.push_back(8 * i + 4);
		indices.push_back(8 * i + 5);

		indices.push_back(8 * i);
		indices.push_back(8 * i + 1);
		indices.push_back(8 * i + 3);

		indices.push_back(8 * i);
		indices.push_back(8 * i + 2);
		indices.push_back(8 * i + 3);

		indices.push_back(8 * i);
		indices.push_back(8 * i + 4);
		indices.push_back(8 * i + 6);

		indices.push_back(8 * i);
		indices.push_back(8 * i + 2);
		indices.push_back(8 * i + 6);

		indices.push_back(8 * i + 7);
		indices.push_back(8 * i + 5);
		indices.push_back(8 * i + 1);

		indices.push_back(8 * i + 7);
		indices.push_back(8 * i + 3);
		indices.push_back(8 * i + 1);

		indices.push_back(8 * i + 7);
		indices.push_back(8 * i + 5);
		indices.push_back(8 * i + 4);

		indices.push_back(8 * i + 7);
		indices.push_back(8 * i + 6);
		indices.push_back(8 * i + 4);

		indices.push_back(8 * i + 7);
		indices.push_back(8 * i + 3);
		indices.push_back(8 * i + 2);

		indices.push_back(8 * i + 7);
		indices.push_back(8 * i + 6);
		indices.push_back(8 * i + 2);
	}
	/*
	indices.push_back((unsigned int)(grid_.N_fluid_ * 8));
	indices.push_back((unsigned int)(grid_.N_fluid_ * 8 + 1));
	indices.push_back((unsigned int)(grid_.N_fluid_ * 8 + 2));
	indices.push_back((unsigned int)(grid_.N_fluid_ * 8 + 3));
	indices.push_back((unsigned int)(grid_.N_fluid_ * 8 + 1));
	indices.push_back((unsigned int)(grid_.N_fluid_ * 8 + 2));
	*/
	return indices;
}

int Fluid_Euler::show(float dt, int n) {
	return draw<Fluid_Euler>(this, dt, n);
}

