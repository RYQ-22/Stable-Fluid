#include "solver.h"


Solver::Solver(int n1, int n2, int n3) :ni(n1), nj(n2), nk(n3) {
	Adiag.resize(n1, n2, n3);
	Aplusi.resize(n1, n2, n3);
	Aplusj.resize(n1, n2, n3);
	Aplusk.resize(n1, n2, n3);
	precon.resize(n1, n2, n3);
	q.resize(n1, n2, n3);
	z.resize(n1, n2, n3);
	s.resize(n1, n2, n3);
}


template <typename T>
T Solver::dot(const Field3<T>& v1, const Field3<T>& v2) {
	if (v1.N1() != v2.N1() || v1.N2() != v2.N2() || v1.N3() != v2.N3()) {
		std::cerr << "2 Vec3 have diffrent size!" << std::endl;
	}
	T ans = 0;
	int size = v1.N1() * v1.N2() * v1.N3();
	for (int idx = 0; idx < size; idx++) {
		ans += v1[idx] * v2[idx];
	}	
	return ans;
}

template <typename T>
void Solver::applyA(const Field3<T>& x, Field3<T>& ans) {
	for (int i = 0; i < x.N1(); i++) for (int j = 0; j < x.N2(); j++) for (int k = 0; k < x.N3(); k++) {
		ans(i, j, k) = Adiag(i, j, k) * x(i, j, k)
			+ Aplusi(i, j, k) * x(i + 1, j, k)
			+ Aplusj(i, j, k) * x(i, j + 1, k)
			+ Aplusk(i, j, k) * x(i, j, k + 1)
			+ Aplusi(i - 1, j, k) * x(i - 1, j, k)
			+ Aplusj(i, j - 1, k) * x(i, j - 1, k)
			+ Aplusk(i, j, k - 1) * x(i, j, k - 1);
	}
	return;
}

void Solver::set_Adiag(float a_diag, const int& i, const int& j, const int& k) {
	Adiag(i, j, k) = a_diag;
}

void Solver::set_Aplusi(float a_plusi, const int& i, const int& j, const int& k) {
	Aplusi(i, j, k) = a_plusi;
}

void Solver::set_Aplusj(float a_plusj, const int& i, const int& j, const int& k) {
	Aplusj(i, j, k) = a_plusj;
}

void Solver::set_Aplusk(float a_plusk, const int& i, const int& j, const int& k) {
	Aplusk(i, j, k) = a_plusk;
}

void Solver::set_zero() {
	Adiag.assign(0);
	Aplusi.assign(0);
	Aplusj.assign(0);
	Aplusk.assign(0);	
}


void Solver::applyPrecon(const Field3<float>& d, const Field3<int>& valid) {
	// compute precon
	float tau = 0.97f;
	float rho = 0.25f;
	float e = 0.0f;
	for (int i = 0; i < ni; i++) for (int j = 0; j < nj; j++) for (int k = 0; k < nk; k++) {
		if (valid(i, j, k) == 1) {
			e = Adiag(i, j, k);
			e += -pow(Aplusi(i - 1, j, k) * precon(i - 1, j, k), 2) - tau * Aplusi(i - 1, j, k) * (Aplusj(i - 1, j, k) + Aplusk(i - 1, j, k)) * pow(precon(i - 1, j, k), 2);
			e += -pow(Aplusj(i, j - 1, k) * precon(i, j - 1, k), 2) - tau * Aplusj(i, j - 1, k) * (Aplusi(i, j - 1, k) + Aplusk(i, j - 1, k)) * pow(precon(i, j - 1, k), 2);
			e += -pow(Aplusk(i, j, k - 1) * precon(i, j, k - 1), 2) - tau * Aplusk(i, j, k - 1) * (Aplusi(i, j, k - 1) + Aplusj(i, j, k - 1)) * pow(precon(i, j, k - 1), 2);
			if (e < rho * Adiag(i, j, k))
				e = Adiag(i, j, k);
			precon(i, j, k) = 1 / sqrt(e + 1e-10f);
		}
		else {
			precon(i, j, k) = 0;
		}
	}
	// solve Lq = d
	float t = 0.0f;
	for (int i = 0; i < ni; i++) for (int j = 0; j < nj; j++) for (int k = 0; k < nk; k++) {
		if (valid(i, j, k)) {
			t = d(i, j, k);
			if (valid(i - 1, j, k))
				t += -Aplusi(i - 1, j, k) * precon(i - 1, j, k) * q(i - 1, j, k);
			if (valid(i, j - 1, k))
				t += -Aplusj(i, j - 1, k) * precon(i, j - 1, k) * q(i, j - 1, k);
			if (valid(i, j, k - 1))
				t += -Aplusk(i, j, k - 1) * precon(i, j, k - 1) * q(i, j, k - 1);
			q(i, j, k) = t * precon(i, j, k);
		}
	}
	//solve L^T z = q
	for (int i = ni - 1; i >= 0; i--) for (int j = nj - 1; j >= 0; j--) for (int k = nk - 1; k >= 0; k--) {
		if (valid(i, j, k)) {
			t = q(i, j, k);
			if (valid(i + 1, j, k))
				t += -Aplusi(i, j, k) * precon(i, j, k) * z(i + 1, j, k);
			if (valid(i, j + 1, k))
				t += -Aplusj(i, j, k) * precon(i, j, k) * z(i, j + 1, k);
			if (valid(i, j, k + 1))
				t += -Aplusk(i, j, k) * precon(i, j, k) * z(i, j, k + 1);
			z(i, j, k) = t * precon(i, j, k);
		}
	}
}


void Solver::solve(Field3<float>& d, Field3<float>& p, const Field3<int>& valid, int max_iterations) {
	// set p = 0
	for (int i = 0; i < ni; i++) for (int j = 0; j < nj; j++) for (int k = 0; k < nk; k++) {
		p(i, j, k) = 0.0f;
	}
	if (d.max() == 0)
		return;
	// PCG
	applyPrecon(d, valid);
	s = z;
	float sigma = dot(z, d);
	float sigma_new;
	float alpha;
	for (int i = 0; i < max_iterations; i++) {
		applyA(s, z);
		alpha = sigma / dot(z, s);
		p += s * alpha;
		d -= z * alpha;
		if (d.max() <= 1e-6f)
			return;
		applyPrecon(d, valid);
		sigma_new = dot(z, d);
		alpha = sigma_new / sigma;
		sigma = sigma_new;
		s = z + s * alpha;
	}
	return;
}
