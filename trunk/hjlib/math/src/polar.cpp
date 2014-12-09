#include "polar.h"

#include <limits>
#include <string.h>

#include "blas_lapack.h"
#include <zjucad/matrix/blas.h>
using namespace zjucad::matrix;

namespace hj {

template <typename T>
static inline void rot2(T pp, T pq, T qq, T &c, T &s)
{
	T t, tao, tao2;
	tao = (qq-pp)/(2*pq);
	tao2 = sqrt(1+tao*tao);
	if(tao >= 0)
		t = 1/(tao+tao2);
	else
		t = 1/(tao-tao2);
	c = 1/sqrt(1+t*t);
	s = t*c;
}

template <typename T, int N>
class Jacobi_eig
{
public:
	Jacobi_eig() {
		V0_ = eye<T>(N);
		row.resize(N);
		col.resize(N);
	}
	int operator()(matrix<T> &A, matrix<T> &V, int steps, T eps) const {
		memcpy(&V[0], &V0_[0], sizeof(T)*N*N);
		int i, j, k;
		bool large_off_diagnal = true;
		for(k = 0; k < steps && large_off_diagnal; ++k) {
			large_off_diagnal = false;
			for(j = 1; j < N; ++j) {
				for(i = 0; i < j; ++i)
					if(rot(A, V, i, j, eps))
						large_off_diagnal = true;
			}
		}

		return k;
	}
private:
	bool rot(matrix<T> &A, matrix<T> &V, int p, int q, T eps) const {
		assert(p < q);
		// cal coef
		const T pq = A(p, q);
		if(fabs(pq) < eps) return false;

//		assert(fabs(pq-A(q, p)) < 1e-4f);

		const T pp = A(p, p), qq = A(q, q);
		T c, s;
		rot2(pp, pq, qq, c, s);
#if 0
		J_ = eye<T>(A.size(1));
		J_(p, p) = c;
		J_(p, q) = s;
		J_(q, p) = -s;
		J_(q, q) = c;
		A = temp(trans(J_)*A);
		A = temp(A*J_);
		V = temp(V*J_);
#else
		rot_JA(A, p, q, c, -s);
		rot_AJ(A, p, q, c, s);
		rot_AJ(V, p, q, c, s);
#endif
		return true;
	}

	inline void rot_JA(matrix<T> &A, int p, int q, T c, T s) const {
#if 0
		row = A(p, colon())*c + A(q, colon())*s;
		A(q, colon()) = A(q, colon())*c - A(p, colon())*s;
		A(p, colon()) = row;
#else
		T *ppA = &A[0];
		T *pA = ppA+p, *qA = ppA+q;
		T tmp;
		for(int iN = 0; iN < N*N; iN+=N) {
			tmp = pA[iN]*c + qA[iN]*s;
			qA[iN] = qA[iN]*c - pA[iN]*s;
			pA[iN] = tmp;
		}
#endif
	}

	inline void rot_AJ(matrix<T> &A, int p, int q, T c, T s) const {
#if 0
		col = A(colon(), p)*c - A(colon(), q)*s;
		A(colon(), q) = A(colon(), p)*s + A(colon(), q)*c;
		A(colon(), p) = col;
#else
		T *ppA = &A[0];
		T *pA = ppA+p*N, *qA = ppA+q*N;
		T tmp;
		for(int i = 0; i < N; ++i) {
			tmp = pA[i]*c - qA[i]*s;
			qA[i] = pA[i]*s + qA[i]*c;
			pA[i] = tmp;
		}
#endif
	}

	mutable matrix<T> A_, J_, row, col;
	matrix<T> V0_;
};

template <typename T, int N>
struct polar_ctx
{
	polar_ctx() {
		V_.resize(N, N);
	}
	matrix<T> ATA_, V_, H_;
	Jacobi_eig<T, N> eig_;
};

template <typename T>
inline void cross(const T *v0, const T *v1, T *v2) {
	v2[0] = v0[1]*v1[2]-v0[2]*v1[1];
	v2[1] = v0[2]*v1[0]-v0[0]*v1[2];
	v2[2] = v0[0]*v1[1]-v0[1]*v1[0];
}


template <typename T>
static T inline det33(const T *m33)
{
#define M(i, j) m33[i+j*3]
	return M(0, 0)*M(1, 1)*M(2, 2)+M(1, 0)*M(2, 1)*M(0, 2)+M(2, 0)*M(0, 1)*M(1, 2)
		-M(0, 0)*M(2, 1)*M(1, 2)-M(1, 0)*M(0, 1)*M(2, 2)-M(2, 0)*M(1, 1)*M(0, 2);
#undef M
}

template <typename T>
int inline min3(T a, T b, T c)
{
	if(a < b) {
		if(b < c)
			return 0;
		else {
			if(a < c)
				return 0;
			else
				return 2;
		}
	}
	else {
		if(a < c)
			return 1;
		else {
			if(b < c)
				return 1;
			else
				return 2;
		}
	}
}

template <typename T>
static int handle_degeneration(polar_ctx<T, 3> *p)
{
	int zero_vec = -1;
	static const int diag_idx[3] = {0, 4, 8};
	for(int i = 0; i < 3; ++i) {
		T diag = p->ATA_[diag_idx[i]];
		const static T eps = std::numeric_limits<T>::epsilon();
		if(diag <= eps) {
			if(zero_vec >= 0) return -1;	// fail
			zero_vec = i;
			continue; 
		}
		p->H_(colon(), i) /= sqrt(diag);
	}
	if(zero_vec >= 0) {
		const int a = (zero_vec+1)%3, b = (zero_vec+2)%3;
		T *pH = &p->H_[0];
		cross(pH+a*3, pH+b*3, pH+zero_vec*3);
	}
	return (zero_vec < 0)?0:1;
}

template <typename T>
static void inline ensure_positive(polar_ctx<T, 3> *p, const T *A, int H_reflected, int reflect)
{
	if(reflect == 2) {
		if(det33(A) < 0 && !H_reflected)
			reflect = 1;
	}
	if(reflect == 1) {
		int axis = min3(p->ATA_[0], p->ATA_[4], p->ATA_[8]);
		T *pH = &p->H_[0]+axis*3;
		pH[0] = -pH[0];
		pH[1] = -pH[1];
		pH[2] = -pH[2];
	}
}

polar3f::polar3f()
{
	ctx_ = new polar_ctx<float, 3>;
}

polar3f::~polar3f()
{
	polar_ctx<float, 3> *p = static_cast<polar_ctx<float, 3> *>(ctx_);
	delete p;
}

static int inline polar3f0(polar_ctx<float, 3> *p, const matrix<float> &A, int reflect, int steps, float eps)
{
	assert(A.size(1) == 3 && A.size(2) == 3 && eps > 0);
	gemm(true, A, false, A, p->ATA_);
	steps = p->eig_(p->ATA_, p->V_, steps, eps);
	gemm(false, A, false, p->V_, p->H_);

	int H_reflected = handle_degeneration(p);
	if(H_reflected < 0) return -1;
	ensure_positive(p, &A[0], H_reflected, reflect);

	return steps;
}

int polar3f::operator()(matrix<float> &A, int reflect, int steps, float eps) const
{
	polar_ctx<float, 3> *p = static_cast<polar_ctx<float, 3> *>(ctx_);
	steps = polar3f0(p, A, reflect, steps, eps);
	if(steps < 0) return -1;
	gemm(false, p->H_, true, p->V_, A);
	return steps;
}

int  polar3f::operator()(const matrix<float> &A, matrix<float> &R, int reflect,
	int steps, float eps) const
{
	polar_ctx<float, 3> *p = static_cast<polar_ctx<float, 3> *>(ctx_);
	steps = polar3f0(p, A, reflect, steps, eps);
	if(steps < 0) return -1;
	gemm(false, p->H_, true, p->V_, R);
	return steps;
}

////

polar3d::polar3d()
{
	ctx_ = new polar_ctx<double, 3>;
}

polar3d::~polar3d()
{
	polar_ctx<double, 3> *p = static_cast<polar_ctx<double, 3> *>(ctx_);
	delete p;
}

static int inline polar3d0(polar_ctx<double, 3> *p, const matrix<double> &A, int reflect, int steps, double eps)
{
	assert(A.size(1) == 3 && A.size(2) == 3 && eps > 0);
	gemm(true, A, false, A, p->ATA_);
	steps = p->eig_(p->ATA_, p->V_, steps, eps);
	gemm(false, A, false, p->V_, p->H_);

	int H_reflected = handle_degeneration(p);
	if(H_reflected < 0) return -1;
	ensure_positive(p, &A[0], H_reflected, reflect);

	return steps;
}

int polar3d::operator()(matrix<double> &A, int reflect, int steps, double eps) const
{
	polar_ctx<double, 3> *p = static_cast<polar_ctx<double, 3> *>(ctx_);
	steps = polar3d0(p, A, reflect, steps, eps);
	if(steps < 0) return -1;
	gemm(false, p->H_, true, p->V_, A);
	return steps;
}

int  polar3d::operator()(const matrix<double> &A, matrix<double> &R, int reflect, int steps, double eps) const
{
	polar_ctx<double, 3> *p = static_cast<polar_ctx<double, 3> *>(ctx_);
	steps = polar3d0(p, A, reflect, steps, eps);
	if(steps < 0) return -1;
	gemm(false, p->H_, true, p->V_, R);
	return steps;
}

///
template <typename T>
class Jacobi_eig2
{
public:
	int operator()(matrix<T> &A, matrix<T> &S, int steps, T eps) const {
		matrix<T> &V = A;
		A_ = A;
		rot(A_, V, eps);	// NOTE: only ONCE!!!
		S[0] = A_[0]; S[1] = A_[3];

		return 1;
	}
private:
	bool rot(matrix<T> &A, matrix<T> &V, T eps) const {
		// cal coef
		T pq = A[2];
		if(fabs(pq) < eps) return false;

//		assert(fabs(pq-A[1]) < 1e-4f);

		T pp = A[0], qq = A[3];
		T c, s;
		rot2(pp, pq, qq, c, s);

		T spq = s*pq, cpq = c*pq;
		A[0] = c*pp-spq;
		A[1] = s*pp+cpq;
		A[2] = cpq-s*qq;
		A[3] = spq+c*qq;

		A[0] = c*A[0]-s*A[2];
		A[3] = s*A[1]+c*A[3];
		// not need this
		//A[1] = 0;
		//A[2] = 0;

		V[0] = V[3] = c;
		V[1] = -s;	V[2] = s;
		return true;
	}
	mutable matrix<T> A_, J_;
};

template <typename T>
struct polar_ctx2
{
	polar_ctx2():ATA_(2, 2), H_(2, 2) {
		S_.resize(2);
	}
	matrix<T> ATA_, S_, H_;
	Jacobi_eig2<T> eig_;
};

template <typename A, typename B, typename C>
inline void mult22_t_nt(const A &a, const B &b, C &c) {
	c[0] = a[0]*b[0]+a[1]*b[1];
	c[1] = a[2]*b[0]+a[3]*b[1];
	c[2] = a[0]*b[2]+a[1]*b[3];
	c[3] = a[2]*b[2]+a[3]*b[3];
}

template <typename A, typename B, typename C>
inline void mult22_nt_nt(const A &a, const B &b, C &c) {
	c[0] = a[0]*b[0]+a[2]*b[1];
	c[1] = a[1]*b[0]+a[3]*b[1];
	c[2] = a[0]*b[2]+a[2]*b[3];
	c[3] = a[1]*b[2]+a[3]*b[3];
}

template <typename A, typename B, typename C>
inline void mult22_nt_t(const A &a, const B &b, C &c) {
	c[0] = a[0]*b[0]+a[2]*b[2];
	c[1] = a[1]*b[0]+a[3]*b[2];
	c[2] = a[0]*b[1]+a[2]*b[3];
	c[3] = a[1]*b[1]+a[3]*b[3];
}

polar2f::polar2f()
{
	ctx_ = new polar_ctx2<float>;
}

polar2f::~polar2f()
{
	polar_ctx2<float> *p = static_cast<polar_ctx2<float> *>(ctx_);
	delete p;
}

int polar2f::operator()(zjucad::matrix::matrix<float> &A, int steps, float eps) const
{
	assert(A.size(1) == 2 && A.size(2) == 2 && eps > 0);
	polar_ctx2<float> *p = static_cast<polar_ctx2<float> *>(ctx_);
	mult22_t_nt(A, A, p->ATA_);
	steps = p->eig_(p->ATA_, p->S_, steps, eps);
	mult22_nt_nt(A, p->ATA_, p->H_);

	float t = sqrt(p->S_[0]);
	p->H_[0] /= t;	p->H_[1] /= t;
	t = sqrt(p->S_[1]);
	p->H_[2] /= t;	p->H_[3] /= t;

	mult22_nt_t(p->H_, p->ATA_, A);
	return steps;
}

///
polar2d::polar2d()
{
	ctx_ = new polar_ctx2<double>;
}

polar2d::~polar2d()
{
	polar_ctx2<double> *p = static_cast<polar_ctx2<double> *>(ctx_);
	delete p;
}

int polar2d::operator()(zjucad::matrix::matrix<double> &A, int steps, double eps) const
{
	assert(A.size(1) == 2 && A.size(2) == 2 && eps > 0);
	polar_ctx2<double> *p = static_cast<polar_ctx2<double> *>(ctx_);
	mult22_t_nt(A, A, p->ATA_);
	steps = p->eig_(p->ATA_, p->S_, steps, eps);
	mult22_nt_nt(A, p->ATA_, p->H_);

	double t = sqrt(p->S_[0]);
	p->H_[0] /= t;	p->H_[1] /= t;
	t = sqrt(p->S_[1]);
	p->H_[2] /= t;	p->H_[3] /= t;

	mult22_nt_t(p->H_, p->ATA_, A);
	return steps;
}

}
