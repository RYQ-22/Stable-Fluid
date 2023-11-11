#include"fluid-Euler.h"

Grid::Grid() {}

Grid::Grid(int N1, int N2, int N3, float l, std::vector<float> phi, std::vector<int> solid) : N1_(N1), N2_(N2), N3_(N3), l_(l) {
	if (phi.size() != N1 * N2 * N3) {
		std::cout << "init level set error!" << std::endl;
		abort();
	}
	if (solid.size() != N1 * N2 * N3) {
		std::cout << "init solid set error!" << std::endl;
		abort();
	}
	x_.resize(N1, N2, N3);
	y_.resize(N1, N2, N3);
	z_.resize(N1, N2, N3);
	Vx_.resize((N1 + 1), N2, N3);
	Vy_.resize(N1, (N2 + 1), N3);
	Vz_.resize(N1, N2, (N3 + 1));
	idx_.resize(N1, N2, N3);
	Vx_valid_.resize(N1 + 1, N2, N3);
	Vy_valid_.resize(N1, N2 + 1, N3);
	Vz_valid_.resize(N1, N2, N3 + 1);
	Adiag.resize(N1, N2, N3);
	Aplusi.resize(N1, N2, N3);
	Aplusj.resize(N1, N2, N3);
	Aplusk.resize(N1, N2, N3);
	phi_ = Vec3<float>(N1, N2, N3, phi);
	solid_ = Vec3<int>(N1, N2, N3, solid);
	for (int i = 0; i < N1; i++) {
		for (int j = 0; j < N2; j++) {
			for (int k = 0; k < N3; k++) {
				// init position
				x_(i, j, k) = i * l_;
				y_(i, j, k) = j * l_;
				z_(i, j, k) = k * l_;
				// init N_fluid_
				if (solid_(i, j, k) == 0 && phi_(i, j, k) < 0.0f) {
					N_fluid_++;
				}
			}
		}
	}
}

float Grid::Vx_ave(float i, float j, float k) {
	int i0 = i > 0 ? (int)i : (int)i - 1;
	int j0 = j > 0 ? (int)j : (int)j - 1;
	int k0 = k > 0 ? (int)k : (int)k - 1;
	float a1 = i0 + 1.0f - i, a2 = i - i0;
	float b1 = j0 + 1.0f - j, b2 = j - j0;
	float c1 = k0 + 1.0f - k, c2 = k - k0;
	return Vx_(i0, j0, k0) * a1 * b1 * c1
		+ Vx_(i0 + 1, j0 + 1, k0 + 1) * a2 * b2 * c2
		+ Vx_(i0 + 1, j0, k0) * a2 * b1 * c1
		+ Vx_(i0, j0 + 1, k0) * a1 * b2 * c1
		+ Vx_(i0, j0, k0 + 1) * a1 * b1 * c2
		+ Vx_(i0 + 1, j0 + 1, k0) * a2 * b2 * c1
		+ Vx_(i0 + 1, j0, k0 + 1) * a2 * b1 * c2
		+ Vx_(i0, j0 + 1, k0 + 1) * a1 * b2 * c2;
}

float Grid::Vy_ave(float i, float j, float k) {
	int i0 = i > 0 ? (int)i : (int)i - 1;
	int j0 = j > 0 ? (int)j : (int)j - 1;
	int k0 = k > 0 ? (int)k : (int)k - 1;
	float a1 = i0 + 1.0f - i, a2 = i - i0;
	float b1 = j0 + 1.0f - j, b2 = j - j0;
	float c1 = k0 + 1.0f - k, c2 = k - k0;
	return Vy_(i0, j0, k0) * a1 * b1 * c1
		+ Vx_(i0 + 1, j0 + 1, k0 + 1) * a2 * b2 * c2
		+ Vy_(i0 + 1, j0, k0) * a2 * b1 * c1
		+ Vy_(i0, j0 + 1, k0) * a1 * b2 * c1
		+ Vy_(i0, j0, k0 + 1) * a1 * b1 * c2
		+ Vy_(i0 + 1, j0 + 1, k0) * a2 * b2 * c1
		+ Vy_(i0 + 1, j0, k0 + 1) * a2 * b1 * c2
		+ Vy_(i0, j0 + 1, k0 + 1) * a1 * b2 * c2;
}

float Grid::Vz_ave(float i, float j, float k) {
	int i0 = i > 0 ? (int)i : (int)i - 1;
	int j0 = j > 0 ? (int)j : (int)j - 1;
	int k0 = k > 0 ? (int)k : (int)k - 1;
	float a1 = i0 + 1.0f - i, a2 = i - i0;
	float b1 = j0 + 1.0f - j, b2 = j - j0;
	float c1 = k0 + 1.0f - k, c2 = k - k0;
	return Vz_(i0, j0, k0) * a1 * b1 * c1
		+ Vz_(i0 + 1, j0 + 1, k0 + 1) * a2 * b2 * c2
		+ Vz_(i0 + 1, j0, k0) * a2 * b1 * c1
		+ Vz_(i0, j0 + 1, k0) * a1 * b2 * c1
		+ Vz_(i0, j0, k0 + 1) * a1 * b1 * c2
		+ Vz_(i0 + 1, j0 + 1, k0) * a2 * b2 * c1
		+ Vz_(i0 + 1, j0, k0 + 1) * a2 * b1 * c2
		+ Vz_(i0, j0 + 1, k0 + 1) * a1 * b2 * c2;
}

float Grid::phi_ave(float i, float j, float k) {
	int i0 = i > 0 ? (int)i : (int)i - 1;
	int j0 = j > 0 ? (int)j : (int)j - 1;
	int k0 = k > 0 ? (int)k : (int)k - 1;
	float a1 = i0 + 1.0f - i, a2 = i - i0;
	float b1 = j0 + 1.0f - j, b2 = j - j0;
	float c1 = k0 + 1.0f - k, c2 = k - k0;
	return phi_(i0, j0, k0) * a1 * b1 * c1
		+ phi_(i0 + 1, j0 + 1, k0 + 1) * a2 * b2 * c2
		+ phi_(i0 + 1, j0, k0) * a2 * b1 * c1
		+ phi_(i0, j0 + 1, k0) * a1 * b2 * c1
		+ phi_(i0, j0, k0 + 1) * a1 * b1 * c2
		+ phi_(i0 + 1, j0 + 1, k0) * a2 * b2 * c1
		+ phi_(i0 + 1, j0, k0 + 1) * a2 * b1 * c2
		+ phi_(i0, j0 + 1, k0 + 1) * a1 * b2 * c2;
}

Fluid_Euler::Fluid_Euler(int N1, int N2, int N3, float l, std::vector<float> phi, std::vector<int> solid) : l_(l) {
	grid_ = Grid(N1, N2, N3, l, phi, solid);
}

void Fluid_Euler::add_force(float dt) {
	for (int i = 0; i < grid_.N1_ + 1; i++) for (int j = 0; j < grid_.N2_ + 1; j++) for (int k = 0; k < grid_.N3_ + 1; k++) {
		// Vy
		if (i < grid_.N1_ && k < grid_.N3_) {
			if ((j < grid_.N2_ && (grid_.phi_(i, j, k) < 0.0f && grid_.solid_(i, j, k) == 0)) || (j != 0 && (grid_.phi_(i, j - 1, k) < 0.0f && grid_.solid_(i, j - 1, k) == 0))) {
				grid_.Vy_(i, j, k) += -G_ * dt;
			}
		}
	}
	return;
}

void Fluid_Euler::advect(float dt) {
	Vec3<float> Vx_temp = grid_.Vx_;
	Vec3<float> Vy_temp = grid_.Vy_;
	Vec3<float> Vz_temp = grid_.Vz_;

	for (int i = 0; i < grid_.N1_ + 1; i++) {
		for (int j = 0; j < grid_.N2_ + 1; j++) {
			for (int k = 0; k < grid_.N3_ + 1; k++) {
				// Vx
				if (j < grid_.N2_ && k < grid_.N3_) {
					if ((i < grid_.N1_ && (grid_.phi_(i, j, k) < 0.0f && grid_.solid_(i, j, k) == 0)) || (i != 0 && (grid_.phi_(i - 1, j, k) < 0.0f && grid_.solid_(i - 1, j, k) == 0))) {
						Vx_temp(i, j, k) = grid_.Vx_ave((float)i - grid_.Vx_(i, j, k) * dt, (float)j - grid_.Vy_ave(i - 0.5f, j + 0.5f, k) * dt, (float)k - grid_.Vz_ave(i - 0.5f, j, k + 0.5f) * dt);
					}
				}
				// Vy
				if (i < grid_.N1_ && k < grid_.N3_) {
					if ((j < grid_.N2_ && (grid_.phi_(i, j, k) < 0.0f && grid_.solid_(i, j, k) == 0)) || (j != 0 && (grid_.phi_(i, j - 1, k) < 0.0f && grid_.solid_(i, j - 1, k) == 0))) {
						Vy_temp(i, j, k) = grid_.Vy_ave((float)i - grid_.Vx_ave(i + 0.5f, j - 0.5f, k) * dt, (float)j - grid_.Vy_(i, j, k) * dt, (float)k - grid_.Vz_ave(i, j - 0.5f, k + 0.5f) * dt);
					}
				}
				// Vz
				if (i < grid_.N1_ && j < grid_.N2_) {
					if ((k < grid_.N3_ && (grid_.phi_(i, j, k) < 0.0f && grid_.solid_(i, j, k) == 0)) || (k != 0 && (grid_.phi_(i, j, k - 1) < 0.0f && grid_.solid_(i, j, k - 1) == 0))) {
						Vz_temp(i, j, k) = grid_.Vx_ave((float)i - grid_.Vx_ave(i + 0.5f, j, k - 0.5f) * dt, (float)j - grid_.Vy_ave(i, j + 0.5f, k - 0.5f) * dt, (float)k - grid_.Vz_(i, j, k) * dt);
					}
				}
			}
		}
	}
	grid_.Vx_ = Vx_temp;
	grid_.Vy_ = Vy_temp;
	grid_.Vz_ = Vz_temp;
	return;
}

void Fluid_Euler::project(float dt) {
	// compute number of fluid grids
	grid_.N_fluid_ = 0;
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
		if (grid_.solid_(i, j, k) == 0 && grid_.phi_(i, j, k) < 0.0f) {
			grid_.idx_(i, j, k) = grid_.N_fluid_;
			grid_.N_fluid_++;
		}
	}

	Eigen::VectorXf p(grid_.N_fluid_);
	Eigen::VectorXf d(grid_.N_fluid_);
	Eigen::SparseMatrix<float> A(grid_.N_fluid_, grid_.N_fluid_);

	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
		if (grid_.solid_(i, j, k) == 0 && grid_.phi_(i, j, k) < 0.0f) {
			int neighbor = 0;
			float di = 0.0f;
			int index = grid_.idx_(i, j, k);
			for (int i1 = -1; i1 <= 1; i1 += 2) {
				neighbor += 1 - grid_.solid_(i + i1, j, k) + 1 - grid_.solid_(i, j + i1, k) + 1 - grid_.solid_(i, j, k + i1);
				if (grid_.solid_(i + i1, j, k) == 0 && grid_.phi_(i + i1, j, k) < 0.0f) {// fluid neighbor
					A.insert(index, grid_.idx_(i + i1, j, k)) = -1.0f;
				}
				if (grid_.solid_(i, j + i1, k) == 0 && grid_.phi_(i, j + i1, k) < 0.0f) {
					A.insert(index, grid_.idx_(i, j + i1, k)) = -1.0f;
				}
				if (grid_.solid_(i, j, k + i1) == 0 && grid_.phi_(i, j, k + i1) < 0.0f) {
					A.insert(index, grid_.idx_(i, j, k + i1)) = -1.0f;
				}
				// calculate divergence of velocity
				if (i1 == 1) {
					di += grid_.solid_(i + 1, j, k) == 0 ? -grid_.Vx_(i + 1, j, k) : 0;
					di += grid_.solid_(i, j + 1, k) == 0 ? -grid_.Vy_(i, j + 1, k) : 0;
					di += grid_.solid_(i, j, k + 1) == 0 ? -grid_.Vz_(i, j, k + 1) : 0;
				}
				if (i1 == -1) {
					di += grid_.solid_(i - 1, j, k) == 0 ? grid_.Vx_(i, j, k) : 0;
					di += grid_.solid_(i, j - 1, k) == 0 ? grid_.Vy_(i, j, k) : 0;
					di += grid_.solid_(i, j, k - 1) == 0 ? grid_.Vz_(i, j, k) : 0;
				}
			}
			A.insert(index, index) = neighbor;
			d(index) = di;
		}
	}

	A.makeCompressed();
	d = d * (1.0f * l_ / dt * rho_);

	// CG
	Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> cg;
	//cg.setMaxIterations(10);
	cg.compute(A);
	p = cg.solve(d);
	// update Vx, Vy, Vz
	//std::cout << grid_.Vy_(10,10,10) << std::endl;

	for (int i = 0; i < grid_.N1_ + 1; i++) {
		for (int j = 0; j < grid_.N2_ + 1; j++) {
			for (int k = 0; k < grid_.N3_ + 1; k++) {
				// Vx
				if (j < grid_.N2_ && k < grid_.N3_) {
					if ((i < grid_.N1_ && (grid_.phi_(i, j, k) < 0.0f && grid_.solid_(i, j, k) == 0)) || (i != 0 && (grid_.phi_(i - 1, j, k) < 0.0f && grid_.solid_(i - 1, j, k) == 0))) {
						if (!(i < grid_.N1_ && grid_.solid_(i, j, k) == 1) && !(i != 0 && grid_.solid_(i - 1, j, k) == 1)) {// not solid
							if (i < grid_.N1_ && grid_.phi_(i, j, k) < 0.0f)// p+
								grid_.Vx_(i, j, k) += -p(grid_.idx_(i, j, k)) / l_ * dt / rho_;
							if (i != 0 && grid_.phi_(i - 1, j, k) < 0.0f)// p-
								grid_.Vx_(i, j, k) += p(grid_.idx_(i - 1, j, k)) / l_ * dt / rho_;
						}
					}
				}
				// Vy
				if (i < grid_.N1_ && k < grid_.N3_) {
					if ((j < grid_.N2_ && (grid_.phi_(i, j, k) < 0.0f && grid_.solid_(i, j, k) == 0)) || (j != 0 && (grid_.phi_(i, j - 1, k) < 0.0f && grid_.solid_(i, j - 1, k) == 0))) {
						if (!(j < grid_.N2_ && grid_.solid_(i, j, k) == 1) && !(j != 0 && grid_.solid_(i, j - 1, k) == 1)) {// not solid
							if (j < grid_.N2_ && grid_.phi_(i, j, k) < 0.0f)// p+
								grid_.Vy_(i, j, k) += -p(grid_.idx_(i, j, k)) / l_ * dt / rho_;
							if (j != 0 && grid_.phi_(i, j - 1, k) < 0.0f)// p-
								grid_.Vy_(i, j, k) += p(grid_.idx_(i, j - 1, k)) / l_ * dt / rho_;
						}
					}
				}
				// Vz
				if (i < grid_.N1_ && j < grid_.N2_) {
					if ((k < grid_.N3_ && (grid_.phi_(i, j, k) < 0.0f && grid_.solid_(i, j, k) == 0)) || (k != 0 && (grid_.phi_(i, j, k - 1) < 0.0f && grid_.solid_(i, j, k - 1) == 0))) {
						if (!(k < grid_.N3_ && grid_.solid_(i, j, k) == 1) && !(k != 0 && grid_.solid_(i, j, k - 1) == 1)) {// not solid
							if (k < grid_.N3_ && grid_.phi_(i, j, k) < 0.0f)// p+
								grid_.Vz_(i, j, k) += -p(grid_.idx_(i, j, k)) / l_ * dt / rho_;
							if (k != 0 && grid_.phi_(i, j, k - 1) < 0.0f)// p-
								grid_.Vz_(i, j, k) += p(grid_.idx_(i, j, k - 1)) / l_ * dt / rho_;
						}
					}
				}
			}
		}
	}
}

void Fluid_Euler::compute_phi(float dt) {
	Vec3<float> phi_temp;
	phi_temp.resize(grid_.N1_, grid_.N2_, grid_.N3_);
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
		phi_temp(i, j, k) = grid_.phi_ave((float)i - grid_.Vx_ave(i + 0.5f, j, k) * dt, (float)j - grid_.Vy_ave(i, j + 0.5f, k) * dt, (float)k - grid_.Vz_ave(i, j, k + 0.5f) * dt);
	}
	grid_.phi_ = phi_temp;

	/*------------------------------------------------- Question ----------------------------------------------------------------
	//extend phi into solid
	for (int i = 0; i < grid_.N1_; i++) for (int j = 0; j < grid_.N2_; j++) for (int k = 0; k < grid_.N3_; k++) {
		if (grid_.phi_(i, j, k) < 0.5f * l_ && grid_.solid_(i,j,k) == 1) {
			grid_.phi_(i, j, k) = -0.5f * l_;
		}
	}
	----------------------------------------------------------------------------------------------------------------------------*/
}

void Fluid_Euler::extrapolate() {
	// reset valid_
	for (int i = 0; i < grid_.N1_ + 1; i++) for (int j = 0; j < grid_.N2_ + 1; j++) for (int k = 0; k < grid_.N3_ + 1; k++) {
		if (j < grid_.N2_ && k < grid_.N3_) {
			if ((i < grid_.N1_ && grid_.solid_(i, j, k) == 0 && grid_.phi_(i, j, k) < 0.0f) || (i - 1 >= 0 && grid_.solid_(i - 1, j, k) == 0 && grid_.phi_(i - 1, j, k) < 0.0f)) {
				grid_.Vx_valid_(i, j, k) = 1;
			}
			else {
				grid_.Vx_valid_(i, j, k) = 0;
			}
		}
		if (i < grid_.N1_ && k < grid_.N3_) {
			if ((j < grid_.N2_ && grid_.solid_(i, j, k) == 0 && grid_.phi_(i, j, k) < 0.0f) || (j - 1 >= 0 && grid_.solid_(i, j - 1, k) == 0 && grid_.phi_(i, j - 1, k) < 0.0f)) {
				grid_.Vy_valid_(i, j, k) = 1;
			}
			else {
				grid_.Vy_valid_(i, j, k) = 0;
			}
		}
		if (i < grid_.N1_ && j < grid_.N2_) {
			if ((k < grid_.N3_ && grid_.solid_(i, j, k) == 0 && grid_.phi_(i, j, k) < 0.0f) || (k - 1 >= 0 && grid_.solid_(i, j, k - 1) == 0 && grid_.phi_(i, j, k - 1) < 0.0f)) {
				grid_.Vz_valid_(i, j, k) = 1;
			}
			else {
				grid_.Vz_valid_(i, j, k) = 0;
			}
		}
	}
	// Apply several iterations of a very simple propagation of valid velocity data in all directions
	for (int num = 0; num < 10; num++) {
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
		if (grid_.solid_(i, j, k) == 0 && grid_.phi_(i, j, k) < 0.0f) {// for liquid
			if ((i + 1 < grid_.N1_ && grid_.solid_(i + 1, j, k) == 1 && grid_.Vx_(i + 1, j, k) > 0))
				grid_.Vx_(i + 1, j, k) = 0;
			if (i > 0 && grid_.solid_(i - 1, j, k) == 1 && grid_.Vx_(i, j, k) < 0)
				grid_.Vx_(i, j, k) = 0;
			if ((j + 1 < grid_.N2_ && grid_.solid_(i, j + 1, k) == 1 && grid_.Vy_(i, j + 1, k) > 0))
				grid_.Vy_(i, j + 1, k) = 0;
			if (j > 0 && grid_.solid_(i, j - 1, k) == 1 && grid_.Vy_(i, j, k) < 0)
				grid_.Vy_(i, j, k) = 0;
			if ((k + 1 < grid_.N3_ && grid_.solid_(i, j, k) == 1 && grid_.Vz_(i, j, k + 1) > 0))
				grid_.Vz_(i, j, k + 1) = 0;
			if (k > 0 && grid_.solid_(i, j, k - 1) == 1 && grid_.Vz_(i, j, k) < 0)
				grid_.Vz_(i, j, k) = 0;
		}
	}
}

void Fluid_Euler::update(float dt) {
	// 1. extrapolation
	extrapolate();
	// 2. update phi
	compute_phi(dt);
	// 3. update v
	// advect
	advect(dt);
	// add force
	add_force(dt);
	// project A * p = d
	project(dt);
	// constrain for fluid
	constrain();// cause unsymmetry
}


std::vector<float> Fluid_Euler::vertices() {
	std::vector<float> vertices;

	grid_.N_fluid_ = 0;
	for (int i0 = 0; i0 < grid_.N1_; i0++) for (int j0 = 0; j0 < grid_.N2_; j0++) for (int k0 = 0; k0 < grid_.N3_; k0++) {
		if (grid_.solid_(i0, j0, k0) == 0 && grid_.phi_(i0, j0, k0) < 0.0f) {
			grid_.N_fluid_++;
			for (int i = -1; i <= 1; i += 2) {
				for (int j = -1; j <= 1; j += 2) {
					for (int k = -1; k <= 1; k += 2) {
						vertices.push_back(10 * (grid_.x_(i0, j0, k0) + 0.5f * l_ * i));
						vertices.push_back(10 * (grid_.y_(i0, j0, k0) + 0.5f * l_ * j));
						vertices.push_back(10 * (grid_.z_(i0, j0, k0) + 0.5f * l_ * k));
						vertices.push_back(0.67843f);
						vertices.push_back(0.84706f);
						vertices.push_back(0.90196f);
					}
				}
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

