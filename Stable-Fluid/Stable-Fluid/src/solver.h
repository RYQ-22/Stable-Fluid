#ifndef SOLVER_H
#define SOLVER_H

#include "vector3.h"
#include "util.h"

template <typename T>
T dot(const Vec3<T>& v1, const Vec3<T>& v2) {
	if (v1.N1() != v2.N1() || v1.N2() != v2.N2() || v1.N3() != v2.N3()) {
		std::cerr << "2 Vec3 have diffrent size!" << std::endl;
	}
	T ans = 0;
	int n1 = v1.N1();
	int n2 = v1.N2();
	int n3 = v1.N3();
	for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++) {		
		ans += v1(i, j, k) * v2(i, j, k);
	}
	return ans;
}

template <typename T>
void applyA(const Vec3<T>& x, Vec3<T>& ans) {
	for (int i = 0; i < x.N1(); i++) for (int j = 0; j < x.N2(); j++) for (int k = 0; k < x.N3(); k++) {
		ans(i, j, k) = grid_.Adiag(i, j, k) * x(i, j, k)
			+ grid_.Aplusi(i, j, k) * x(i + 1, j, k)
			+ grid_.Aplusj(i, j, k) * x(i, j + 1, k)
			+ grid_.Aplusk(i, j, k) * x(i, j, k + 1)
			+ grid_.Aplusi(i - 1, j, k) * x(i - 1, j, k)
			+ grid_.Aplusj(i, j - 1, k) * x(i, j - 1, k)
			+ grid_.Aplusk(i, j, k - 1) * x(i, j, k - 1);			
	}
	return;
}

class Solver { // solve Ap = d
public:
	Solver(int n1, int n2, int n3);
	void applyPrecon();
	void solve(Vec3<float>& d, Vec3<float>& p);
	void set_Adiag(int);
	void set_Aplusi(int);
	void set_Aplusj(int);
	void set_Aplusk(int);

private:
	int ni, nj, nk;
	Vec3<int> Adiag, Aplusi, Aplusj, Aplusk;
	Vec3<float> d, p;
	// z: auxiliary vetor
	// s: search vector
	Vec3<float> precon, q, z, s;
};

Solver::Solver(int n1, int n2, int n3) :ni(n1), nj(n2), nk(n3) {
	Adiag.resize(n1, n2, n3);
	Aplusi.resize(n1, n2, n3);
	Aplusj.resize(n1, n2, n3);
	Aplusk.resize(n1, n2, n3);
	d.resize(n1, n2, n3);
	p.resize(n1, n2, n3);
	precon.resize(n1, n2, n3);
	q.resize(n1, n2, n3);
	z.resize(n1, n2, n3);
	s.resize(n1, n2, n3);
}

/*
void Solver::applyPrecon() {
	// compute precon
	float tau = 0.97f;
	float rho = 0.25f;
	float e = 0.0f;
	for (int i = 0; i < ni; i++) for (int j = 0; j < nj; j++) for (int k = 0; k < nk; k++) {
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
*/
#endif
