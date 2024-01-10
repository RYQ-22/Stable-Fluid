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
	new_Vx_.resize((N1 + 1), N2, N3);
	new_Vy_.resize(N1, (N2 + 1), N3);
	new_Vz_.resize(N1, N2, (N3 + 1));
	new_Vx_aux_.resize((N1 + 1), N2, N3);
	new_Vy_aux_.resize(N1, (N2 + 1), N3);
	new_Vz_aux_.resize(N1, N2, (N3 + 1));
	Vx_reflected_.resize((N1 + 1), N2, N3);
	Vy_reflected_.resize(N1, (N2 + 1), N3);
	Vz_reflected_.resize(N1, N2, (N3 + 1));

	idx_.resize(N1, N2, N3);

	Vx_valid_.resize(N1 + 1, N2, N3);
	Vy_valid_.resize(N1, N2 + 1, N3);
	Vz_valid_.resize(N1, N2, N3 + 1);
	old_Vx_valid_.resize(N1 + 1, N2, N3);
	old_Vy_valid_.resize(N1, N2 + 1, N3);
	old_Vz_valid_.resize(N1, N2, N3 + 1);

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

	new_phi_.resize(N1, N2, N3);
	phi_ = Field3<float>(N1, N2, N3, phi);
	solid_phi_ = Field3<float>(N1, N2, N3, solid_phi);

	rigidBody_phi_.resize(N1, N2, N3);
	new_solid_phi_.resize(N1, N2, N3);

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

float Grid::solid_phi_ave(const glm::vec3& pos) {
	return interpolate_value(pos.x, pos.y, pos.z, solid_phi_);
}

float Grid::rigidBody_phi_ave(float i, float j, float k) {
	return interpolate_value(i, j, k, rigidBody_phi_);
}

float Grid::rigidBody_phi_ave(const glm::vec3& pos) {
	return interpolate_value(pos.x, pos.y, pos.z, rigidBody_phi_);
}

Fluid_Euler::Fluid_Euler(int N1, int N2, int N3, float l, std::vector<float> phi, std::vector<float> solid_phi) : l_(l) {
	grid_ = Grid(N1, N2, N3, l, phi, solid_phi);

	// make the particles large enough so they always appear on the grid
	particle_radius_ = (float)(l_ * 1.01 * sqrt(3.0) / 2.0);
	// init particles
	int seed = 0;
	for (int n = 0; n < 64*4; n++) {
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

void Fluid_Euler::set_velocity(glm::vec3 v) {
	grid_.Vx_.assign(v.x);
	grid_.Vy_.assign(v.y);
	grid_.Vz_.assign(v.z);
	return;
}

// Apply RK2 to advect a point in the domain.
glm::vec3 Fluid_Euler::trace_rk2(const glm::vec3& position, float dt) {
	glm::vec3 input = position;
	glm::vec3 velocity = get_velocity(input);
	velocity = get_velocity(input + 0.5f * dt * velocity);
	input += dt * velocity;
	return input;
}

void Fluid_Euler::semi_Lagrangian(const Field3<float>& field, Field3<float>& new_field, float dt) {	
	#pragma omp parallel for schedule(dynamic)
	{
		for (int i = 0; i < field.N1(); i++) for (int j = 0; j < field.N2(); j++) for (int k = 0; k < field.N3(); k++) {
			glm::vec3 pos = glm::vec3(i * l_, j * l_, k * l_);
			new_field(i, j, k) = interpolate_value(trace_rk2(pos, -dt) / l_, field);
		}
	}
	
	return;
}

void Fluid_Euler::mac_Cormack(const Field3<float>& field, Field3<float>& new_field, Field3<float>& new_field_aux, float dt) {		
	// BFECC
	semi_Lagrangian(field, new_field, dt);
	semi_Lagrangian(new_field, new_field_aux, -dt);
	new_field = new_field + (field - new_field_aux) * 0.5f;
	glm::vec3 pos, source_pos;
	float min_val, max_val;
	for (int i = 0; i < field.N1(); i++) for (int j = 0; j < field.N2(); j++) for (int k = 0; k < field.N3(); k++) {
		pos = glm::vec3(i * l_, j * l_, k * l_);
		source_pos = trace_rk2(pos, -dt);
		min_val = sample_min(source_pos / l_, field);
		max_val = sample_max(source_pos / l_, field);
		if (new_field(i, j, k) < min_val || new_field(i, j, k) > max_val) {
			new_field(i, j, k) = interpolate_value(source_pos, field);
		}
	}	
}

glm::vec3 Fluid_Euler::get_velocity(const glm::vec3& position) {
	float vx_ave = grid_.Vx_ave(position.x / l_ + 0.5f, position.y / l_, position.z / l_);
	float vy_ave = grid_.Vy_ave(position.x / l_, position.y / l_ + 0.5f, position.z / l_);
	float vz_ave = grid_.Vz_ave(position.x / l_, position.y / l_, position.z / l_ + 0.5f);
	return glm::vec3(vx_ave, vy_ave, vz_ave);
}

bool Fluid_Euler::valid(int i, int j, int k) {// fluid, and not in solid
	if (add_rigidBody_) {
		if (i >= 0 && i < grid_.N1_ && j >= 0 && j < grid_.N2_ && k >= 0 && k < grid_.N3_ && grid_.solid_phi_(i, j, k) > 0 && grid_.phi_(i, j, k) < 0.0f && grid_.rigidBody_phi_(i, j, k) > 0) {
			return true;
		}
	}
	else {
		if (i >= 0 && i < grid_.N1_ && j >= 0 && j < grid_.N2_ && k >= 0 && k < grid_.N3_ && grid_.solid_phi_(i, j, k) > 0 && grid_.phi_(i, j, k) < 0.0f) {
			return true;
		}
	}
	return false;
}

template <typename T>
T Fluid_Euler::dot(const Field3<T>& v1, const Field3<T>& v2) {
	if (v1.N1() != v2.N1() || v1.N2() != v2.N2() || v1.N3() != v2.N3()) {
		std::cerr << "2 Vec3 have diffrent size!" << std::endl;
	}
	T ans = 0;
	int N = v1.size();
	#pragma omp parallel for
	for (int idx = 0; idx < N; idx++) {
		ans += v1[idx] * v2[idx];
	}
	return ans;
}

void Fluid_Euler::applyA(const Field3<float>& x, Field3<float>& ans) {	
	#pragma omp parallel for
	for (int i = 0; i < x.N1(); i++) for (int j = 0; j < x.N2(); j++) for (int k = 0; k < x.N3(); k++) {
		ans(i, j, k) = grid_.Adiag(i, j, k) * x(i, j, k);		
		ans(i, j, k) += grid_.Aplusi(i, j, k) * x(i + 1, j, k);
		ans(i, j, k) += grid_.Aplusj(i, j, k) * x(i, j + 1, k);
		ans(i, j, k) += grid_.Aplusk(i, j, k) * x(i, j, k + 1);
		ans(i, j, k) += grid_.Aplusi(i - 1, j, k) * x(i - 1, j, k);
		ans(i, j, k) += grid_.Aplusj(i, j - 1, k) * x(i, j - 1, k);
		ans(i, j, k) += grid_.Aplusk(i, j, k - 1) * x(i, j, k - 1);
	}
	if (add_rigidBody_&0) {
		float a = rho_ * l_ * l_ * l_ / rigidBody_.M();
		glm::mat3 b = rho_ * l_ * l_ * l_ * rigidBody_.I_inv();
		grid_.apply_J.assign(0);
		grid_.apply_J += grid_.J_x * dot(grid_.J_x, x) * a;
		grid_.apply_J += grid_.J_y * dot(grid_.J_y, x) * a;
		grid_.apply_J += grid_.J_z * dot(grid_.J_z, x) * a;
		grid_.apply_J += grid_.Jrot_x * (dot(grid_.Jrot_x, x) * b[0][0] + dot(grid_.Jrot_y, x) * b[0][1] + dot(grid_.Jrot_z, x) * b[0][2]);
		grid_.apply_J += grid_.Jrot_y * (dot(grid_.Jrot_x, x) * b[1][0] + dot(grid_.Jrot_y, x) * b[1][1] + dot(grid_.Jrot_z, x) * b[1][2]);
		grid_.apply_J += grid_.Jrot_z * (dot(grid_.Jrot_x, x) * b[2][0] + dot(grid_.Jrot_y, x) * b[2][1] + dot(grid_.Jrot_z, x) * b[2][2]);
		ans -= grid_.apply_J;
	}
	return;
}

Field3<float> Fluid_Euler::applyA(const Field3<float>& x) {
	Field3<float> ans(x.N1(), x.N2(), x.N3());
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
	if (add_rigidBody_) {
		float a = rho_ * l_ * l_ * l_ / rigidBody_.M();
		glm::mat3 b = rho_ * l_ * l_ * l_ * rigidBody_.I_inv();
		grid_.apply_J.assign(0);
		grid_.apply_J += grid_.J_x * dot(grid_.J_x, x) * a;
		grid_.apply_J += grid_.J_y * dot(grid_.J_y, x) * a;
		grid_.apply_J += grid_.J_z * dot(grid_.J_z, x) * a;
		grid_.apply_J += grid_.Jrot_x * (dot(grid_.Jrot_x, x) * b[0][0] + dot(grid_.Jrot_y, x) * b[0][1] + dot(grid_.Jrot_z, x) * b[0][2]);
		grid_.apply_J += grid_.Jrot_y * (dot(grid_.Jrot_x, x) * b[1][0] + dot(grid_.Jrot_y, x) * b[1][1] + dot(grid_.Jrot_z, x) * b[1][2]);
		grid_.apply_J += grid_.Jrot_z * (dot(grid_.Jrot_x, x) * b[2][0] + dot(grid_.Jrot_y, x) * b[2][1] + dot(grid_.Jrot_z, x) * b[2][2]);
		ans -= grid_.apply_J;
	}
	return ans;
}

void Fluid_Euler::add_force(float dt) {
	// Vy	
	#pragma omp parallel
	{
		for (int idx = 0; idx < grid_.Vy_.size(); idx++) {
			grid_.Vy_[idx] += -G_ * dt;
			if (use_reflection_)
				grid_.new_Vy_[idx] += -G_ * dt;
		}
	}
	return;
}


void Fluid_Euler::advect_particles(float dt) {
	//glm::vec3 pos_1, pos_2;	
	for (int i = 0; i < particles_.size(); i++) {
		// BFECC
		//pos_1 = trace_rk2(particles_[i], dt);
		//pos_2 = trace_rk2(pos_1, -dt);
		//particles_[i] = pos_1 + (particles_[i] - pos_2) / 2.0f;
		particles_[i] = trace_rk2(particles_[i], dt);
		// boundary solid
		float phi_val = grid_.solid_phi_ave(particles_[i] / l_);
		if (phi_val < 0) {
			glm::vec3 grad;
			interpolate_gradient(grad, particles_[i] / l_, grid_.solid_phi_);
			grad = glm::normalize(grad);
			particles_[i] -= phi_val * grad;
		}
		// boundary rigidBody
		if (add_rigidBody_) {
			float phi_val = grid_.rigidBody_phi_ave(particles_[i] / l_);
			if (phi_val < 0) {
				glm::vec3 grad;
				interpolate_gradient(grad, particles_[i] / l_, grid_.rigidBody_phi_);
				grad = glm::normalize(grad);
				particles_[i] -= phi_val * grad;
			}
		}		
	}
}


void Fluid_Euler::advect(float dt, int field_id) {
	if (use_mc_) {
		// BFECC
		if (field_id == 0) {
			mac_Cormack(grid_.Vx_, grid_.new_Vx_, grid_.new_Vx_aux_, dt);
			mac_Cormack(grid_.Vy_, grid_.new_Vy_, grid_.new_Vy_aux_, dt);
			mac_Cormack(grid_.Vz_, grid_.new_Vz_, grid_.new_Vz_aux_, dt);
		}
		else if (field_id == 1) {
			mac_Cormack(grid_.Vx_reflected_, grid_.new_Vx_, grid_.new_Vx_aux_, dt);
			mac_Cormack(grid_.Vy_reflected_, grid_.new_Vy_, grid_.new_Vy_aux_, dt);
			mac_Cormack(grid_.Vz_reflected_, grid_.new_Vz_, grid_.new_Vz_aux_, dt);
		}
		else {
			abort();
		}
							
	}
	else {
		// semi-Lagrangian
		if (field_id == 0) {
			semi_Lagrangian(grid_.Vx_, grid_.new_Vx_, dt);
			semi_Lagrangian(grid_.Vy_, grid_.new_Vy_, dt);
			semi_Lagrangian(grid_.Vz_, grid_.new_Vz_, dt);
		}
		else if (field_id == 1) {
			semi_Lagrangian(grid_.Vx_reflected_, grid_.new_Vx_, dt);
			semi_Lagrangian(grid_.Vy_reflected_, grid_.new_Vy_, dt);
			semi_Lagrangian(grid_.Vz_reflected_, grid_.new_Vz_, dt);
		}
		else {
			abort();
		}
	}	
	grid_.Vx_ = grid_.new_Vx_;
	grid_.Vy_ = grid_.new_Vy_;
	grid_.Vz_ = grid_.new_Vz_;	
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
	grid_.d.assign(0);
	#pragma omp parallel for
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {		
		if (valid(i, j, k)) {// for liquid
			grid_.Adiag(i, j, k) = 0;
			grid_.Aplusi(i, j, k) = 0;
			grid_.Aplusj(i, j, k) = 0;
			grid_.Aplusk(i, j, k) = 0;			
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
		else {
			grid_.Adiag(i, j, k) = 0;
			grid_.Aplusi(i, j, k) = 0;
			grid_.Aplusj(i, j, k) = 0;
			grid_.Aplusk(i, j, k) = 0;
		}
	}
	// solve: A * p = d
	solve(15);
	
	if (add_surfaceTension_) {
		for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_ + 1; j++) for (int k = 0; k < grid_.N3_ + 1; k++) {
			if (grid_.solid_phi_(i, j, k) > 0 && grid_.phi_(i, j, k) > 0) {
				grid_.p(i, j, k) = surfaceTensorCoefficient_ / l_ / l_ / l_ * 0.001 / rho_ *
					(-6 * grid_.phi_(i, j, k) + grid_.phi_(i - 1, j, k) + grid_.phi_(i + 1, j, k) + grid_.phi_(i, j - 1, k) + grid_.phi_(i, j + 1, k) + grid_.phi_(i, j, k - 1) + grid_.phi_(i, j, k + 1));
			}
		}
	}

	// update Vx, Vy, Vz
	#pragma omp parallel for
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

void Fluid_Euler::apply_pressure_liquid(Field3<float>& vx, Field3<float>& vy, Field3<float>& vz) {
	#pragma omp parallel for
	for (int i = 0; i < grid_.N1_ + 1; i++) for (int j = 0; j < grid_.N2_ + 1; j++) for (int k = 0; k < grid_.N3_ + 1; k++) {
		// Vx
		if (j < grid_.N2_ && k < grid_.N3_) {
			if ((i < grid_.N1_ && valid(i, j, k)) || (i != 0 && valid(i - 1, j, k))) {// liquid neighbor
				if (!(i < grid_.N1_ && grid_.solid_phi_(i, j, k) < 0) && !(i != 0 && grid_.solid_phi_(i - 1, j, k) < 0)) {// not solid
					if (i < grid_.N1_ && grid_.phi_(i, j, k) < 0.0f)// p+
						vx(i, j, k) += -grid_.p(i, j, k);
					if (i != 0 && grid_.phi_(i - 1, j, k) < 0.0f)// p-
						vx(i, j, k) += grid_.p(i - 1, j, k);
				}
			}
		}
		// Vy
		if (i < grid_.N1_ && k < grid_.N3_) {
			if ((j < grid_.N2_ && valid(i, j, k)) || (j != 0 && valid(i, j - 1, k))) {
				if (!(j < grid_.N2_ && grid_.solid_phi_(i, j, k) < 0) && !(j != 0 && grid_.solid_phi_(i, j - 1, k) < 0)) {// not solid
					if (j < grid_.N2_ && grid_.phi_(i, j, k) < 0.0f)// p+
						vy(i, j, k) += -grid_.p(i, j, k);
					if (j != 0 && grid_.phi_(i, j - 1, k) < 0.0f)// p-
						vy(i, j, k) += grid_.p(i, j - 1, k);
				}
			}
		}
		// Vz
		if (i < grid_.N1_ && j < grid_.N2_) {
			if ((k < grid_.N3_ && valid(i, j, k)) || (k != 0 && valid(i, j, k - 1))) {
				if (!(k < grid_.N3_ && grid_.solid_phi_(i, j, k) < 0) && !(k != 0 && grid_.solid_phi_(i, j, k - 1) < 0)) {// not solid
					if (k < grid_.N3_ && grid_.phi_(i, j, k) < 0.0f)// p+
						vz(i, j, k) += -grid_.p(i, j, k);
					if (k != 0 && grid_.phi_(i, j, k - 1) < 0.0f)// p-
						vz(i, j, k) += grid_.p(i, j, k - 1);
				}
			}
		}
	}
}


void Fluid_Euler::project_new() {
	// compute A, d
	grid_.d.assign(0);
	#pragma omp parallel for
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
		/*grid_.Adiag(i, j, k) = computeVolume(i, j, k, grid_.phi_, 1) + computeVolume(i, j, k, grid_.phi_, 2) + computeVolume(i, j, k, grid_.phi_, 3)
			+ computeVolume(i - 1, j, k, grid_.phi_, 1) + computeVolume(i, j - 1, k, grid_.phi_, 2) + computeVolume(i, j, k - 1, grid_.phi_, 3);
		grid_.Aplusi(i, j, k) = -computeVolume(i, j, k, grid_.phi_, 1);
		grid_.Aplusj(i, j, k) = -computeVolume(i, j, k, grid_.phi_, 2);
		grid_.Aplusk(i, j, k) = -computeVolume(i, j, k, grid_.phi_, 3);
		grid_.d(i, j, k) = -computeVolume(i, j, k, grid_.phi_, 1) * grid_.Vx_(i + 1, j, k) + computeVolume(i - 1, j, k, grid_.phi_, 1) * grid_.Vx_(i, j, k)
			- computeVolume(i, j, k, grid_.phi_, 2) * grid_.Vx_(i, j + 1, k) + computeVolume(i, j - 1, k, grid_.phi_, 2) * grid_.Vy_(i, j, k)
			- computeVolume(i, j, k, grid_.phi_, 3) * grid_.Vx_(i, j, k + 1) + computeVolume(i, j, k - 1, grid_.phi_, 3) * grid_.Vz_(i, j, k);
		*/

		if (valid(i, j, k)) {// for liquid
			grid_.Adiag(i, j, k) = 0;
			grid_.Aplusi(i, j, k) = 0;
			grid_.Aplusj(i, j, k) = 0;
			grid_.Aplusk(i, j, k) = 0;			
			for (int i1 = -1; i1 <= 1; i1 += 2) {
				if (grid_.solid_phi_(i + i1, j, k) > 0 && grid_.rigidBody_phi_(i + i1, j, k) > 0) {
					grid_.Adiag(i, j, k)++;
				}
				if (grid_.solid_phi_(i, j + i1, k) > 0 && grid_.rigidBody_phi_(i, j + i1, k) > 0) {
					grid_.Adiag(i, j, k)++;
				}
				if (grid_.solid_phi_(i, j, k + i1) > 0 && grid_.rigidBody_phi_(i, j, k + i1) > 0) {
					grid_.Adiag(i, j, k)++;
				}
				if (i1 == 1) {
					if (grid_.solid_phi_(i + 1, j, k) > 0) {
						if (grid_.rigidBody_phi_(i + 1, j, k) < 0) {
							float vx = (rigidBody_.V(glm::vec3(i + 0.5f, j, k) * l_)).x;
							if (grid_.Vx_(i + 1, j, k) * sign(vx) > vx) {
								grid_.d(i, j, k) += -vx;
							}
							else {
								grid_.d(i, j, k) += -grid_.Vx_(i + 1, j, k);
							}
						}
						else {
							grid_.d(i, j, k) += -grid_.Vx_(i + 1, j, k);
						}
					}
					if (grid_.solid_phi_(i, j + 1, k) > 0) {
						if (grid_.rigidBody_phi_(i, j + 1, k) < 0) {
							float vy = (rigidBody_.V(glm::vec3(i, j + 0.5f, k) * l_)).y;
							if (grid_.Vy_(i, j + 1, k) * sign(vy) > vy) {
								grid_.d(i, j, k) += -vy;
							}
							else {
								grid_.d(i, j, k) += -grid_.Vy_(i, j + 1, k);
							}
						}
						else {
							grid_.d(i, j, k) += -grid_.Vy_(i, j + 1, k);
						}
					}
					if (grid_.solid_phi_(i, j, k + 1) > 0) {
						if (grid_.rigidBody_phi_(i, j, k + 1) < 0) {
							float vz = (rigidBody_.V(glm::vec3(i, j, k + 0.5f) * l_)).z;
							if (grid_.Vz_(i, j, k + 1) * sign(vz) > vz) {
								grid_.d(i, j, k) += -vz;
							}
							else {
								grid_.d(i, j, k) += -grid_.Vz_(i, j, k + 1);
							}
						}
						else {
							grid_.d(i, j, k) += -grid_.Vz_(i, j, k + 1);
						}
					}					
				}
				if (i1 == -1) {
					if (grid_.solid_phi_(i - 1, j, k) > 0) {
						if (grid_.rigidBody_phi_(i - 1, j, k) < 0) {
							float vx = (rigidBody_.V(glm::vec3(i - 0.5f, j, k) * l_)).x;
							if (grid_.Vx_(i, j, k) * sign(vx) < vx) {
								grid_.d(i, j, k) += vx;
							}
							else {
								grid_.d(i, j, k) += grid_.Vx_(i, j, k);
							}
						}
						else {
							grid_.d(i, j, k) += grid_.Vx_(i, j, k);
						}
					}
					if (grid_.solid_phi_(i, j - 1, k) > 0) {
						if (grid_.rigidBody_phi_(i, j - 1, k) < 0) {
							float vy = (rigidBody_.V(glm::vec3(i, j - 0.5f, k) * l_)).y;
							if (grid_.Vy_(i, j, k) * sign(vy) < vy) {
								grid_.d(i, j, k) += vy;

							}
							else {
								grid_.d(i, j, k) += grid_.Vy_(i, j, k);
							}
						}
						else {
							grid_.d(i, j, k) += grid_.Vy_(i, j, k);
						}
					}
					if (grid_.solid_phi_(i, j, k - 1) > 0) {
						if (grid_.rigidBody_phi_(i, j, k - 1) < 0) {
							float vz = (rigidBody_.V(glm::vec3(i, j, k - 0.5f) * l_)).z;
							if (grid_.Vz_(i, j, k) * sign(vz) < vz) {
								grid_.d(i, j, k) += vz;
							}
							else {
								grid_.d(i, j, k) += grid_.Vz_(i, j, k);
							}
						}
						else {
							grid_.d(i, j, k) += grid_.Vz_(i, j, k);
						}
					}
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
		else {
			grid_.Adiag(i, j, k) = 0;
			grid_.Aplusi(i, j, k) = 0;
			grid_.Aplusj(i, j, k) = 0;
			grid_.Aplusk(i, j, k) = 0;
		}
	}
	
	// compute J
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
		grid_.J_x(i, j, k) = computeVolume(i, j, k, grid_.rigidBody_phi_, 1) - computeVolume(i - 1, j, k, grid_.rigidBody_phi_, 1);
		grid_.J_y(i, j, k) = computeVolume(i, j, k, grid_.rigidBody_phi_, 2) - computeVolume(i, j - 1, k, grid_.rigidBody_phi_, 2);
		grid_.J_z(i, j, k) = computeVolume(i, j, k, grid_.rigidBody_phi_, 3) - computeVolume(i, j, k - 1, grid_.rigidBody_phi_, 3);
		grid_.Jrot_x(i, j, k) = (computeVolume(i, j, k, grid_.rigidBody_phi_, 1) - computeVolume(i, j, k - 1, grid_.rigidBody_phi_, 1)) * (j * l_ - rigidBody_.c().y)
			- (computeVolume(i, j, k, grid_.rigidBody_phi_, 1) - computeVolume(i, j - 1, k, grid_.rigidBody_phi_, 1)) * (k * l_ - rigidBody_.c().z);
		grid_.Jrot_y(i, j, k) = (computeVolume(i, j, k, grid_.rigidBody_phi_, 2) - computeVolume(i - 1, j, k, grid_.rigidBody_phi_, 2)) * (k * l_ - rigidBody_.c().z)
			- (computeVolume(i, j, k, grid_.rigidBody_phi_, 2) - computeVolume(i, j, k - 1, grid_.rigidBody_phi_, 2)) * (i * l_ - rigidBody_.c().x);
		grid_.Jrot_z(i, j, k) = (computeVolume(i, j, k, grid_.rigidBody_phi_, 3) - computeVolume(i, j - 1, k, grid_.rigidBody_phi_, 3)) * (i * l_ - rigidBody_.c().x)
			- (computeVolume(i, j, k, grid_.rigidBody_phi_, 3) - computeVolume(i - 1, j, k, grid_.rigidBody_phi_, 3)) * (j * l_ - rigidBody_.c().y);
	}

	// d* = d - J^T V / dx^2
	//grid_.d += grid_.J_x * rigidBody_.Vc().x + grid_.J_y * rigidBody_.Vc().y + grid_.J_z * rigidBody_.Vc().z
	//	+ grid_.Jrot_x * rigidBody_.omega().x + grid_.Jrot_y * rigidBody_.omega().y + grid_.Jrot_z * rigidBody_.omega().z;
	

	// solve: A * p = d
	solve(7);

	// update Vx, Vy, Vz
	#pragma omp parallel for
	for (int i = 0; i < grid_.N1_ + 1; i++) for (int j = 0; j < grid_.N2_ + 1; j++) for (int k = 0; k < grid_.N3_ + 1; k++) {
		// Vx
		if (j < grid_.N2_ && k < grid_.N3_) {
			if ((i < grid_.N1_ && valid(i, j, k)) || (i != 0 && valid(i - 1, j, k))) {// liquid neighbor
				if (!(i < grid_.N1_ && grid_.solid_phi_(i, j, k) < 0) && !(i != 0 && grid_.solid_phi_(i - 1, j, k) < 0)) {// not solid
					if (i < grid_.N1_ )// p+
						grid_.Vx_(i, j, k) += -grid_.p(i, j, k);
					if (i != 0 )// p-
						grid_.Vx_(i, j, k) += grid_.p(i - 1, j, k);
				}
			}
		}
		// Vy
		if (i < grid_.N1_ && k < grid_.N3_) {
			if ((j < grid_.N2_ && valid(i, j, k)) || (j != 0 && valid(i, j - 1, k))) {
				if (!(j < grid_.N2_ && grid_.solid_phi_(i, j, k) < 0) && !(j != 0 && grid_.solid_phi_(i, j - 1, k) < 0)) {// not solid
					if (j < grid_.N2_ )// p+
						grid_.Vy_(i, j, k) += -grid_.p(i, j, k);
					if (j != 0 )// p-
						grid_.Vy_(i, j, k) += grid_.p(i, j - 1, k);
				}
			}
		}
		// Vz
		if (i < grid_.N1_ && j < grid_.N2_) {
			if ((k < grid_.N3_ && valid(i, j, k)) || (k != 0 && valid(i, j, k - 1))) {
				if (!(k < grid_.N3_ && grid_.solid_phi_(i, j, k) < 0) && !(k != 0 && grid_.solid_phi_(i, j, k - 1) < 0)) {// not solid
					if (k < grid_.N3_ )// p+
						grid_.Vz_(i, j, k) += -grid_.p(i, j, k);
					if (k != 0 )// p-
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
	// solve L^T z = q
	// Mz = r
	//grid_.z = grid_.d;
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
	grid_.p.assign(0);
	if (grid_.d.max() == 0) {		
		return;
	}
	// PCG
	//applyPrecon();
	grid_.z = grid_.d;

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
		//applyPrecon();		
		grid_.z = grid_.d;

		sigma_new = dot(grid_.z, grid_.d);
		alpha = sigma_new / sigma;
		sigma = sigma_new;
		grid_.s = grid_.z + grid_.s * alpha;
	}
	
	return;
}

void Fluid_Euler::compute_phi(float dt) {	
	grid_.phi_.assign(3 * l_);
	// Vec3<float> phi_temp = grid_.phi_;
	#pragma omp parallel
	{
		for (int p = 0; p < particles_.size(); p++) {
			glm::ivec3 cell_ind(particles_[p] / l_);
			for (int k = max(0, cell_ind[2] - 1); k <= min(cell_ind[2] + 1, grid_.N3_ - 1); ++k) {
				for (int j = max(0, cell_ind[1] - 1); j <= min(cell_ind[1] + 1, grid_.N2_ - 1); ++j) {
					for (int i = max(0, cell_ind[0] - 1); i <= min(cell_ind[0] + 1, grid_.N1_ - 1); ++i) {
						glm::vec3 sample_pos((i)*l_, (j)*l_, (k)*l_);
						float test_val = glm::length(sample_pos - particles_[p]) - particle_radius_;
						if (test_val < grid_.phi_(i, j, k)) {
							grid_.phi_(i, j, k) = test_val;
						}
					}
				}
			}
		}
	}	
	
	// smooth the implicit surface	
	for (int num = 0; num < 1; num++) {
		glm::vec3 phi_grad;
		for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
			interpolate_gradient(phi_grad, (float)i, (float)j, (float)k, grid_.phi_);
			grid_.new_phi_(i, j, k) = -grid_.phi_(i, j, k) / std::sqrt(pow(grid_.phi_(i, j, k), 2) + pow(l_, 2)) * (glm::length(phi_grad) - 1);
		}
		grid_.phi_ = grid_.new_phi_;
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
	extrapolate(grid_.Vx_, grid_.new_Vx_, grid_.Vx_valid_, grid_.old_Vx_valid_);
	extrapolate(grid_.Vy_, grid_.new_Vy_, grid_.Vy_valid_, grid_.old_Vy_valid_);
	extrapolate(grid_.Vz_, grid_.new_Vz_, grid_.Vz_valid_, grid_.old_Vz_valid_);
}

void Fluid_Euler::extrapolate(Field3<float>& v, Field3<float> v_new, Field3<int>& valid, Field3<int>& valid_old) {
	v_new = v;
	for (int num = 0; num < 10; num++) {
		float sum = 0;
		int count = 0;
		valid_old = valid;
		for (int i = 0; i < v.N1(); i++) for (int j = 0; j < v.N2(); j++) for (int k = 0; k < v.N3(); k++) {
			sum = 0;
			count = 0;
			if (!valid_old(i, j, k)) {

				if (valid_old(i + 1, j, k)) {
					sum += v(i + 1, j, k);
					count++;
				}
				if (valid_old(i - 1, j, k)) {
					sum += v(i - 1, j, k);
					count++;
				}
				if (valid_old(i, j + 1, k)) {
					sum += v(i, j + 1, k);
					count++;
				}
				if (valid_old(i, j - 1, k)) {
					sum += v(i, j - 1, k);
					count++;
				}
				if (valid_old(i, j, k + 1)) {
					sum += v(i, j, k + 1);
					count++;
				}
				if (valid_old(i, j, k - 1)) {
					sum += v(i, j, k - 1);
					count++;
				}

				if (count > 0) {
					v_new(i, j, k) = sum / count;
					valid(i, j, k) = 1;
				}
			}
		}
		v = v_new;
	}
}

void Fluid_Euler::constrain() {
	#pragma omp parallel for
	{
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
	if (add_rigidBody_) {
		#pragma omp parallel for
		{
			for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
				if (valid(i, j, k)) {// for liquid
					if (i + 1 < grid_.N1_ && grid_.rigidBody_phi_(i + 1, j, k) < 0) {
						float vx = (rigidBody_.V(glm::vec3(i + 0.5f, j, k) * l_)).x;
						if (grid_.Vx_(i + 1, j, k) * sign(vx) > vx)	grid_.Vx_(i + 1, j, k) = vx;
					}
					if (i > 0 && grid_.rigidBody_phi_(i - 1, j, k) < 0) {
						float vx = (rigidBody_.V(glm::vec3(i - 0.5f, j, k) * l_)).x;
						if (grid_.Vx_(i, j, k) * sign(vx) < vx) grid_.Vx_(i, j, k) = vx;
					}
					if (j + 1 < grid_.N2_ && grid_.rigidBody_phi_(i, j + 1, k) < 0) {
						float vy = (rigidBody_.V(glm::vec3(i, j + 0.5f, k) * l_)).y;
						if (grid_.Vy_(i, j + 1, k) * sign(vy) > vy)	grid_.Vy_(i, j + 1, k) = vy;
					}						
					if (j > 0 && grid_.rigidBody_phi_(i, j - 1, k) < 0) {
						float vy = (rigidBody_.V(glm::vec3(i, j - 0.5f, k) * l_)).y;
						if (grid_.Vy_(i, j, k) * sign(vy) < vy)	grid_.Vy_(i, j, k) = vy;
					}						
					if (k + 1 < grid_.N3_ && grid_.rigidBody_phi_(i, j, k) < 0) {
						float vz = (rigidBody_.V(glm::vec3(i, j, k + 0.5) * l_)).z;
						if (grid_.Vz_(i, j, k + 1) * sign(vz) > vz)	grid_.Vz_(i, j, k + 1) = vz;
					}						
					if (k > 0 && grid_.rigidBody_phi_(i, j, k - 1) < 0) {
						float vz = (rigidBody_.V(glm::vec3(i, j, k - 0.5) * l_)).z;
						if (grid_.Vz_(i, j, k) * sign(vz) < vz)	grid_.Vz_(i, j, k) = vz;
					}						
				}
			}
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
	// debug
	std::cout << "V_liquid: " << max(grid_.Vx_.max(), grid_.Vy_.max(), grid_.Vz_.max()) << " c_z: " << rigidBody_.c().z << " c_y: " << rigidBody_.c().y  << " c_x: " << rigidBody_.c().x << " p_max: " << grid_.p.max() << " d_max: " << grid_.d.max() << std::endl;
	//std::cout << " J_x max: " << grid_.J_x.max() << " J_y max: " << grid_.J_y.max() << " J_z max: " << grid_.J_z.max() << std::endl;
	//std::cout << "Jrot_x max: " << grid_.Jrot_x.max() << " Jrot_y max: " << grid_.Jrot_y.max() << " Jrot_z max: " << grid_.Jrot_z.max() << std::endl;
	//std::cout << grid_.Vy_reflected_(20,20,20) << std::endl;
	
	// 1. extrapolation
	extrapolate();
	// 2. update phi
	advect_particles(dt);
	compute_phi(dt);
	// update rigid body
	if (add_rigidBody_) {
		update_rigidBody(dt);
	}
	// 3. update v
	if (use_reflection_) {
		// advect-reflect solver:
		// u1 = A(u0,u0,dt/2)
		// u2 = P(u1)
		// u3 = 2u2-u1
		// u4 = A(u3,u2,dt/2)
		// u5 = P(u4)
		advect(dt / 2);// new_v_ = A(v_)
		add_force(dt);// new_v_ = A(v_) + f dt/2
		project();// v_ = P(new_v_)
		
		grid_.Vx_reflected_ = grid_.Vx_;
		grid_.Vy_reflected_ = grid_.Vy_;
		grid_.Vz_reflected_ = grid_.Vz_;
		apply_pressure_liquid(grid_.Vx_reflected_, grid_.Vy_reflected_, grid_.Vz_reflected_);
		//grid_.Vx_reflected_ = grid_.Vx_ * 2.0f - grid_.new_Vx_;
		//grid_.Vy_reflected_ = grid_.Vy_ * 2.0f - grid_.new_Vy_;
		//grid_.Vz_reflected_ = grid_.Vz_ * 2.0f - grid_.new_Vz_;
		advect(dt / 2, 1);		
		project();
	}
	else {
		// advect-project solver:
		advect(dt);
		add_force(dt);
		if (!add_rigidBody_) {
			project();
		}
		else {
			project_new();
			if (!use_fixedOperation_) apply_pressure_rigidBody();
		}
	}	
	// constrain for fluid
	constrain();// cause unsymmetry	
}

void Fluid_Euler::add_rigidBody(
	std::string obj_path,
	float scale, 
	glm::vec3 translate,
	float M,
	glm::mat3 I_ref,
	glm::vec3 c0,
	glm::quat rotationQuaternion,
	glm::vec3 Vc,
	glm::vec3 omega) {
	if (add_rigidBody_) {
		std::cerr << "Rigid body has been added." << std::endl;
		abort();
	}
	add_rigidBody_ = true;
	// read from .obj	
	// init rigidBody_phi_
	std::vector<float> phi(grid_.N1_ * grid_.N2_ * grid_.N3_);
	obj_2_SDF_py(grid_.N1_, grid_.N2_, grid_.N3_, l_, obj_path, phi, scale, translate.x, translate.y, translate.z);
	for (int idx = 0; idx < phi.size(); idx++) {
		grid_.rigidBody_phi_[idx] = phi[idx];
	}
	// init rigidBody_
	std::vector<glm::vec3> rigidBody_vertices;
	std::vector<Face> rigidBody_faces;
	read_from_OBJ(obj_path, rigidBody_vertices, rigidBody_faces, scale, translate);
	rigidBody_ = RigidBody(rigidBody_vertices, rigidBody_faces, M, I_ref, c0, rotationQuaternion, Vc, omega);

	grid_.J_x.resize(grid_.N1_, grid_.N2_, grid_.N3_);
	grid_.J_y.resize(grid_.N1_, grid_.N2_, grid_.N3_);
	grid_.J_z.resize(grid_.N1_, grid_.N2_, grid_.N3_);
	grid_.Jrot_x.resize(grid_.N1_, grid_.N2_, grid_.N3_);
	grid_.Jrot_y.resize(grid_.N1_, grid_.N2_, grid_.N3_);
	grid_.Jrot_z.resize(grid_.N1_, grid_.N2_, grid_.N3_);
	grid_.apply_J.resize(grid_.N1_, grid_.N2_, grid_.N3_);
	// init rigidBody_particles
	int seed = 0;
	for (int n = 0; n < 64*4; n++) {
		for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
			float a = randhashf(seed++); float b = randhashf(seed++); float c = randhashf(seed++);
			float x = (float)i + a, y = (float)j + b, z = (float)k + c;
			if (grid_.rigidBody_phi_ave(x, y, z) <= -particle_radius_) {
				if (grid_.solid_phi_ave(x, y, z) > 0)
					rigidBody_particles_.push_back(glm::vec3(x * l_, y * l_, z * l_));
			}
		}
	}
}

void Fluid_Euler::update_rigidBody(float dt) {
	// update rigidBody_
	if (!use_fixedOperation_) {
		rigidBody_.update(dt, grid_.solid_phi_, l_);
	}
	else {
		if (rigidBody_.c().z <= 0.6f) {
			rigidBody_.movePosition(glm::vec3(0.0f, 0.0f, 0.5f) * dt);
			rigidBody_.setVelocity(glm::vec3(0.0f, 0.0f, 0.5f));
		}
	}
	// debug
	// std::cout << rigidBody_.Vc().y << " " << rigidBody_.c().y << std::endl;
	// update rigidBody_phi_
	grid_.rigidBody_phi_.assign(3 * l_);
	#pragma omp parallel
	{
		for (int p = 0; p < rigidBody_particles_.size(); p++) {
			glm::ivec3 cell_ind((rigidBody_.c() + rigidBody_.rotationMatrix() * (rigidBody_particles_[p] - rigidBody_.c0())) / l_);
			for (int k = max(0, cell_ind[2] - 1); k <= min(cell_ind[2] + 1, grid_.N3_ - 1); ++k) {
				for (int j = max(0, cell_ind[1] - 1); j <= min(cell_ind[1] + 1, grid_.N2_ - 1); ++j) {
					for (int i = max(0, cell_ind[0] - 1); i <= min(cell_ind[0] + 1, grid_.N1_ - 1); ++i) {
						glm::vec3 sample_pos((i)*l_, (j)*l_, (k)*l_);
						float test_val = glm::length(sample_pos - (rigidBody_.c() + rigidBody_.rotationMatrix() * (rigidBody_particles_[p] - rigidBody_.c0()))) - particle_radius_;
						if (test_val < grid_.rigidBody_phi_(i, j, k)) {
							grid_.rigidBody_phi_(i, j, k) = test_val;
						}
					}
				}
			}
		}
	}	
	return;
}

void Fluid_Euler::set_fixedOperation() {
	use_fixedOperation_ = true;
}

float Fluid_Euler::computeVolume(int i, int j, int k, const Field3<float>& phi, int id) {
	/* id = 
	1: +x
	2: +y
	3: +z
	4: -x
	5: -y
	6: -z
	*/
	float phi0, phi1, phi2, phi3;
	if (id == 1) {
		phi0 = interpolate_value(i + 0.5f, j + 0.5f, k + 0.5f, phi);
		phi1 = interpolate_value(i + 0.5f, j - 0.5f, k + 0.5f, phi);
		phi2 = interpolate_value(i + 0.5f, j + 0.5f, k - 0.5f, phi);
		phi3 = interpolate_value(i + 0.5f, j - 0.5f, k - 0.5f, phi);
	}
	else if (id == 2) {
		phi0 = interpolate_value(i + 0.5f, j + 0.5f, k + 0.5f, phi);
		phi1 = interpolate_value(i - 0.5f, j + 0.5f, k + 0.5f, phi);
		phi2 = interpolate_value(i + 0.5f, j + 0.5f, k - 0.5f, phi);
		phi3 = interpolate_value(i - 0.5f, j + 0.5f, k - 0.5f, phi);
	}
	else if (id == 3) {
		phi0 = interpolate_value(i + 0.5f, j + 0.5f, k + 0.5f, phi);
		phi1 = interpolate_value(i - 0.5f, j + 0.5f, k + 0.5f, phi);
		phi2 = interpolate_value(i + 0.5f, j - 0.5f, k + 0.5f, phi);
		phi3 = interpolate_value(i - 0.5f, j - 0.5f, k + 0.5f, phi);
	}
	else if (id == 4) {
		phi0 = interpolate_value(i - 0.5f, j + 0.5f, k + 0.5f, phi);
		phi1 = interpolate_value(i - 0.5f, j - 0.5f, k + 0.5f, phi);
		phi2 = interpolate_value(i - 0.5f, j + 0.5f, k - 0.5f, phi);
		phi3 = interpolate_value(i - 0.5f, j - 0.5f, k - 0.5f, phi);
	}
	else if (id == 5) {
		phi0 = interpolate_value(i + 0.5f, j - 0.5f, k + 0.5f, phi);
		phi1 = interpolate_value(i - 0.5f, j - 0.5f, k + 0.5f, phi);
		phi2 = interpolate_value(i + 0.5f, j - 0.5f, k - 0.5f, phi);
		phi3 = interpolate_value(i - 0.5f, j - 0.5f, k - 0.5f, phi);
	}
	else if (id == 6) {
		phi0 = interpolate_value(i + 0.5f, j + 0.5f, k - 0.5f, phi);
		phi1 = interpolate_value(i - 0.5f, j + 0.5f, k - 0.5f, phi);
		phi2 = interpolate_value(i + 0.5f, j - 0.5f, k - 0.5f, phi);
		phi3 = interpolate_value(i - 0.5f, j - 0.5f, k - 0.5f, phi);
	}
	else {
		std::cerr << "invalid input id" << std::endl;
		abort();
	}	
	// compute s1: 0 1 3
	float s1 = computeTriangleArea(phi0, phi1, phi3);
	// compute s2: 0 2 3
	float s2 = computeTriangleArea(phi0, phi2, phi3);
	return (s1 + s2) / 2;
}

void Fluid_Euler::apply_pressure_rigidBody() {
	glm::vec3 f = glm::vec3(0);
	glm::vec3 tau = glm::vec3(0);
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
		f += -computeVolume(i, j, k, grid_.rigidBody_phi_, 1) * (grid_.p(i + 1, j, k) - grid_.p(i, j, k)) * glm::vec3(1.0f, 0.0f, 0.0f);
		f += -computeVolume(i, j, k, grid_.rigidBody_phi_, 2) * (grid_.p(i, j + 1, k) - grid_.p(i, j, k)) * glm::vec3(0.0f, 1.0f, 0.0f);
		f += -computeVolume(i, j, k, grid_.rigidBody_phi_, 3) * (grid_.p(i, j, k + 1) - grid_.p(i, j, k)) * glm::vec3(0.0f, 0.0f, 1.0f);
		tau += computeVolume(i, j, k, grid_.rigidBody_phi_, 1) *
			((grid_.p(i, j + 1, k) - grid_.p(i, j, k)) * (k * l_ - rigidBody_.c().z)
				- (grid_.p(i, j, k + 1) - grid_.p(i, j, k)) * (j * l_ - rigidBody_.c().y)) * glm::vec3(1.0f, 0.0f, 0.0f);
		tau += computeVolume(i, j, k, grid_.rigidBody_phi_, 2) *
			((grid_.p(i, j, k + 1) - grid_.p(i, j, k)) * (i * l_ - rigidBody_.c().x)
				- (grid_.p(i + 1, j, k) - grid_.p(i, j, k)) * (k * l_ - rigidBody_.c().z)) * glm::vec3(0.0f, 1.0f, 0.0f);
		tau += computeVolume(i, j, k, grid_.rigidBody_phi_, 3) *
			((grid_.p(i + 1, j, k) - grid_.p(i, j, k)) * (j * l_ - rigidBody_.c().y)
				- (grid_.p(i, j + 1, k) - grid_.p(i, j, k)) * (i * l_ - rigidBody_.c().x)) * glm::vec3(0.0f, 0.0f, 1.0f);
	}
	f *= rho_ * l_ * l_ * l_ * 3.68;
	tau *= rho_ * l_ * l_ * l_ * 3.0;
	rigidBody_.applyImpulse(f, 1);
	rigidBody_.applyTorqueImpulse(tau, 1);
	return;
}

std::vector<float> Fluid_Euler::vertices() {
	std::vector<float> vertices;

	// liquid
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
	// particles
	for (int p = 0; p < particles_.size(); p++) {
		for (int i = -1; i <= 1; i += 2) for (int j = -1; j <= 1; j += 2) for (int k = -1; k <= 1; k += 2) {
			vertices.push_back(10 * (particles_[p].x + 0.5f * l_ * i));
			vertices.push_back(10 * (particles_[p].y + 0.5f * l_ * j));
			vertices.push_back(10 * (particles_[p].z + 0.5f * l_ * k));
			if (grid_.solid_phi_ave(particles_[p]/l_) <= 0) {
				vertices.push_back(0.67843f);
				vertices.push_back(0.14706f);
				vertices.push_back(0.10196f);
			}
			else {
				vertices.push_back(0.67843f);
				vertices.push_back(0.84706f);
				vertices.push_back(0.90196f);
			}			
		}
	}
	*/
	

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

void Fluid_Euler::run(float dt, int n) {
	while (1) {		
		// output fps
		auto currentTime_chrono = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration<double>(currentTime_chrono.time_since_epoch());
		double currentTime = duration.count();
		deltaTime_ = currentTime - lastTime_;
		lastTime_ = currentTime;
		deltaTime_tot_ += deltaTime_;
		frameRate_ += 1.0 / deltaTime_;
		count_++;
		// print
		if (deltaTime_tot_ >= 0.5) {
			//std::cout << std::setprecision(4) << frameRate_ / count_ << " FPS\n";
			frameRate_ = 0;
			count_ = 0;
			deltaTime_tot_ = 0.0;
		}
		for (int i = 0; i < n; i++) {
			update(dt);
		}
	}	
}

void Fluid_Euler::outputPLY(int fps, int t, float dt, int n) {
	// delete .ply
	std::string folderPath_ = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/ply";
	try {
		for (const auto& entry : std::filesystem::directory_iterator(folderPath_)) {
			std::filesystem::remove(entry.path());
		}
	}
	catch (const std::filesystem::filesystem_error& err) {
		std::cerr << "Error during cleanup: " << err.what() << std::endl;
	}
	// delete sdf_i.txt before quit
	CleanupHelper cleanup("C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/sdf");
	
	std::string fileName;
	for (int frame_num = 0; frame_num < fps * t; frame_num++) {
		// output liquid SDF
		fileName = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/sdf/sdf_" + std::to_string(frame_num) + ".txt";
		std::ofstream file(fileName);
		if (!file.is_open()) {
			std::cerr << "Can not open the file " << fileName << std::endl;
			return;
		}
		file << grid_.N1_ << " " << grid_.N2_ << " " << grid_.N3_ << " " << l_ << "\n";
		for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
			file << grid_.phi_(i, j, k) << "\n";
		}
		file.close();
		std::cout << "sdf " + std::to_string(frame_num) + "\n";
		
		// output rigidBody .ply		
		if (add_rigidBody_) {
			std::string rigidBody_path = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/ply/rigidBody_" + std::to_string(frame_num) + ".ply";
			write_to_PLY(rigidBody_path, rigidBody_.vertices(), rigidBody_.faces(), 10.0f);
		}		
		
		// update system
		for (int i = 0; i < n; i++) {
			update(dt);
		}
	}
	// output liquid .ply
	std::string command;
	for (int frame_num = 0; frame_num < fps * t; frame_num++) {
		// call python script
		command = "D:/Anaconda/envs/myenv/python.exe C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/sdf2ply.py -frame "
			+ std::to_string(frame_num);
		system(command.c_str());
		std::cout << "ply " + std::to_string(frame_num) + "\n";
	}	
	return;
}

void Fluid_Euler::outputXYZ(int fps, int t, float dt, int n) {
	// delete .xyz
	std::string folderPath_ = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/xyz";
	try {
		for (const auto& entry : std::filesystem::directory_iterator(folderPath_)) {
			std::filesystem::remove(entry.path());
		}
	}
	catch (const std::filesystem::filesystem_error& err) {
		std::cerr << "Error during cleanup: " << err.what() << std::endl;
	}
	// output .xyz
	std::string fileName;
	for (int frame_num = 0; frame_num < fps * t; frame_num++) {
		fileName = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/xyz/xyz_" + std::to_string(frame_num) + ".xyz";
		std::ofstream file(fileName);
		if (!file.is_open()) {
			std::cerr << "Can not open the file " << fileName << std::endl;
			return;
		}		
		for (int p = 0; p < particles_.size(); p++) {
			file << particles_[p].x << " " << particles_[p].y << " " << particles_[p].z << "\n";
		}
		file.close();
		std::cout << "xyz " + std::to_string(frame_num) + "\n";
		for (int i = 0; i < n; i++) {
			update(dt);
		}
	}
	return;
}

void Fluid_Euler::simple_render(std::string output_file, int fps) {
	// delete .png before quit
	CleanupHelper cleanup("C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/png");
	// call python script
	std::string command1 = "D:/Anaconda/envs/myenv/python.exe C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/pv_renderer.py";
	std::string command2 = "D:/Anaconda/envs/myenv/python.exe C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/png2mp4.py -output_file "
		+ output_file
		+ " -fps " + std::to_string(fps);
	// use pyvista to render .ply
	system(command1.c_str());
	// png to mp4
	system(command2.c_str());	
	return;
}

void Fluid_Euler::pbrt_render(std::string output_file, int fps, int t) {
	// TO DO...		
	CleanupHelper cleanup("C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/png");	

	std::string pbrt_folder = "C:/Users/11862/Desktop/vs_code/ACG/pbrt-v4/build/Release";
	std::string command1 = "D:/Anaconda/envs/myenv/python.exe C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/pbrt_renderer.py -pbrt_folder "
		+ pbrt_folder
		+ " -frame " + std::to_string(fps * t);
	if (add_rigidBody_) {
		command1 = command1 + " -add_rigidBody " + std::to_string(true);
	}
	std::string command2 = "D:/Anaconda/envs/myenv/python.exe C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/png2mp4.py -output_file "
		+ output_file
		+ " -fps " + std::to_string(fps);
	// use pbrt to render .ply
	system(command1.c_str());
	// png to mp4
	system(command2.c_str());
	// delete .pbrt
	try {
		std::filesystem::path input_folder = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/scene";
		for (const auto& entry : std::filesystem::directory_iterator(input_folder)) { 
			if (entry.is_regular_file() && entry.path().extension() == ".pbrt") {
				std::filesystem::remove(entry.path());
			}
		}
	}
	catch (const std::filesystem::filesystem_error& err) {
		std::cerr << "Error during cleanup: " << err.what() << std::endl;
	}
	return;
}


