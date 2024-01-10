#ifndef SOLVER_H
#define SOLVER_H

#include "field.h"
#include "util.h"

typedef glm::vec<6, float> vec6;
typedef glm::mat<6, 6, float> mat6;

class Solver { // solve Ap = d
public:
	Solver(int n1, int n2, int n3);

	template <typename T>
	T dot(const Field3<T>& v1, const Field3<T>& v2);
	template <typename T>
	void applyA(const Field3<T>& x, Field3<T>& ans);

	void set_Adiag(float a_diag, const int& i, const int& j, const int& k);
	void set_Aplusi(float a_plusi, const int& i, const int& j, const int& k);
	void set_Aplusj(float a_plusj, const int& i, const int& j, const int& k);
	void set_Aplusk(float a_plusk, const int& i, const int& j, const int& k);
	void set_zero();

	void applyPrecon(const Field3<float>& d, const Field3<int>& valid);
	void solve(Field3<float>& d, Field3<float>& p, const Field3<int>& valid, int max_iterations);	

private:
	int ni, nj, nk;
	Field3<float> Adiag, Aplusi, Aplusj, Aplusk;	
	//std::vector<Field3<vec6>> J;
	//std::vector<mat6> Ms_inverse;	
	// z: auxiliary vetor
	// s: search vector
	Field3<float> precon, q, z, s;
};

#endif
