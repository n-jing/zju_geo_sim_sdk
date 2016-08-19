#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <ctime>
#include <memory>

#include <fstream>

#include <zjucad/matrix/include/matrix.h>
#include <zjucad/matrix/include/io.h>

#include "sparse.h"

using namespace hj::sparse;
using namespace std;
using namespace zjucad::matrix;

#define ENABLE_DEPRECATED 0

/*
#include <boost/numeric/ublas/matrix_sparse.hpp>
using namespace boost::numeric;

#include "sparse_multi_cl.h"

void test_boost()
{
	// map<i*rows+j, value> ...
	ublas::mapped_matrix<double> A0(4, 3);

	// not work?
	ublas::mapped_matrix<double, ublas::column_major,
		ublas::map_array<std::size_t, double> > A1(4, 3);

	//csc or csr
	ublas::compressed_matrix<double> A2(4, 3);

	//The values are stored in 3 parallel array as triples (i, j, value).
	//More than one value for each pair of indices is possible,
	//the real value is the sum of all.
	ublas::coordinate_matrix<double> A(4, 3);
	A2(1, 1) = 5;
	A2(1, 2) = 3;
	A2(2, 0) = 2;
	A2(0, 1) = 1;
}
*/

template <typename T>
void write_32(ofstream &ofs, const matrix<T> &v)
{
	int n = v.size();
	ofs.write((const char *)&n, sizeof(int));
	ofs.write((const char *)&v[0], n*sizeof(T));
}

template <typename T>
void read_32(ifstream &ifs, matrix<T> &v)
{
	int n;
	ifs.read((char *)&n, sizeof(int));
	v.resize(n);
	ifs.read((char *)&v[0], n*sizeof(T));
}

// void test_solver()
// {
// 	{ // test the crash at:		rtn->ptr = cholmod_l_analyze(pA, pc);
// 		// ifstream ifs("crash-sp.zjumat", ifstream::binary);
// 		// if(ifs.fail()) {
// 		// 	cerr << "read fail." << endl;
// 		// 	return;
// 		// }
// 		// csc<double> A;
// 		// read(ifs, A.ptr());
// 		// read(ifs, A.idx());
// 		// read(ifs, A.val());
// 		// cout << "read over." << endl;
// 		// A.rows_ = 27;

// 		csc<double, int> A32;
// 		// A32.ptr() = A.ptr();
// 		// A32.idx() = A.idx();
// 		// A32.val() = A.val();
// 		// {
// 		// 	ofstream ofs("crash-sp-int32.zjumat", ofstream::binary);
// 		// 	write_32(ofs, A32.ptr());
// 		// 	write_32(ofs, A32.idx());
// 		// 	write_32(ofs, A32.val());
// 		// }
// 		{
// 			ifstream ifs("crash-sp-int32.zjumat", ifstream::binary);
// 			read_32(ifs, A32.ptr());
// 			read_32(ifs, A32.idx());
// 			read_32(ifs, A32.val());
// 		}
// 		A32.rows_ = 27;

// //		cout << A32.val();
// //		matrix<double> dA;
// //		convert(A32, dA);
// //		cout << dA;
// 		cout << nnz(A32) << endl;
// 		// auto_ptr<hj::sparse::direct_solver_A> sv;
// 		// sv.reset(hj::sparse::direct_solver_A::create(A32, "cholmod"));
// 		auto_ptr<hj::sparse::solver> sv;
// 		sv.reset(hj::sparse::solver::create(A32, "cholmod"));
// 	}
// 	return;

//    	const int DIM = 500, NRHS = 3;
// 	matrix<double> A = to<int>(pow(rand<double>(DIM, DIM), 10)*1.02),
// 		b = floor(rand<double>(DIM, NRHS)*30),
// 		x = floor(rand<double>(DIM, NRHS)*30);
// 	matrix<double> ATA = trans(A)*A+eye<double>(DIM)*1e-5, c;

// 	csc<double, int> sATA;
// 	convert(ATA, sATA, 1e-5);
// 	cout << "nnz ratio: " << nnz(sATA)/double(ATA.size())*100 << endl;

// 	clock_t beg, elapse;

// 	const char *solver_name[] = {
// 		"umfpack",
// 		"cholmod"
// 	};
// 	const int LOOP = 100;
// 	for(int s = 0; s < 2; ++s) {
// #if ENABLE_DEPRECATED
// 		x(colon()) = 0;
// 		beg = clock();
// 		for(int i = 0; i < LOOP; ++i) {
// 			auto_ptr<hj::sparse::solver> sv;
// 			sv.reset(hj::sparse::solver::create(sATA, solver_name[s]));
// 			if(!sv.get()) {
// 				cout << "create fail." << endl;
// 				return;
// 			}
// 			sv->solve(&b[0], &x[0], NRHS);
// 		}
// 		elapse = clock()-beg;
// 		cout << "deprecated: " << solver_name[s] << " " << elapse << endl;
 
// 		c = -b;
// 		spm_mm(false, sATA, false, x, false ,c);
// 		cout << norm(c) << endl;
// #endif
// 		x(colon()) = 0;
// 		beg = clock();
// 		for(int i = 0; i < LOOP; ++i) {
// 			auto_ptr<hj::sparse::direct_solver_A> sv;
// 			sv.reset(hj::sparse::direct_solver_A::create(sATA, solver_name[s]));
// 			if(!sv.get()) {
// 				cout << "create fail." << endl;
// 				return;
// 			}
// 			if(sv->solve(&b[0], &x[0], NRHS))
// 				cout << "solve fail." << endl;
// 		}
// 		elapse = clock()-beg;
// 		cout << "current   : " << solver_name[s] << " " << elapse << endl;
// 		c = -b;
// 		mm(false, sATA, x, c);
// 		cout << norm(c) << endl;
// 	}
// 	cout << "finish test_solver." << endl;
// }

class timer {
public:
	timer(const char *name)
		:name_(name) {
		t_ = clock();
	}
	~timer() {
		t_ = clock()-t_;
		cout << name_ << ":\t" << t_ << endl;
	}
	string name_;
	clock_t t_;
};

int test_basic_op(void)
{
	const int m = 70, n = 50, k = 10;
	matrix<double> A = floor(pow(rand<double>(m, k), 10)*100), b = floor(rand<double>(k, n)*30);
	csc<double> spA;
	convert(A, spA, 1e-4);
	cout << "nnz ratio: " << nnz(spA)/double(A.size()) << endl;
 	matrix<double> c;

	size_t LOOP;

//	cout << A;
//	cout << convert(spA, A);

	cout << "convert" << endl;
	{
		matrix<double> A1;
		convert(spA, A1);
		cout << norm(A1-A) << endl;
	}

	cout << "mm should be zero" << endl;
	LOOP = 20000;
#if ENABLE_DEPRECATED
    c = -A*b;
	{timer t("spm_mm");
		for(size_t i = 0; i < LOOP; ++i)
			spm_mm(false, spA, false, b, false, c);
	}
	c = -A*b;
	cout << norm(spm_mm(false, spA, false, b, false, c)) << endl;
#endif
	c = -A*b;
	{timer t("    mm");
		for(size_t i = 0; i < LOOP; ++i)
			mm(false, spA, b, c);
	}
	c = -A*b;
	cout << norm(mm(false, spA, b, c)) << endl;
#if ENABLE_DEPRECATED
	cout << "mv should be zero" << endl;
	c = -A*b;
	{timer t("spm_mxv");
		for(size_t i = 0; i < LOOP; ++i)
			spm_mxv(spA, b, c);
	}
	c = -A*b;
	spm_mxv(spA, b, c);
	cout << norm(c(colon(), 0)) << endl;
	c = -trans(A)*b;
	{timer t("spm_mTxv");
		for(size_t i = 0; i < LOOP; ++i)
			spm_mTxv(spA, b, c);
	}
	c = -trans(A)*b;
	spm_mTxv(spA, b, c);
	cout << norm(c(colon(), 0)) << endl;
#endif
	c = -A*b;
	{timer t("    mm");
		for(size_t i = 0; i < LOOP; ++i)
			mv(false, spA, b, c);
	}
	c = -A*b;
	mv(false, spA, b, c);
	cout << norm(c(colon(), 0)) << endl;

	c = -trans(A)*b;
	{timer t("    mv");
		for(size_t i = 0; i < LOOP; ++i)
			mv(true, spA, b, c);
	}
	c = -trans(A)*b;
	mv(true, spA, b, c);
	cout << norm(c(colon(), 0)) << endl;

	cout << "AAT should be zero" << endl;
	matrix<double> AAT0 = A*trans(A), AAT1;
	csc<double> spAAT;
	AAT<std_map>(spA, spAAT);
	LOOP /= 50;
	{timer t("std_map");
		for(size_t i = 0; i < LOOP; ++i)
			AAT<std_map>(spA, spAAT);
	}
	{timer t("sorted_vec");
		for(size_t i = 0; i < LOOP; ++i)
			AAT<map_by_sorted_vector>(spA, spAAT);
	}
	{timer t("unsorted_vec");
		for(size_t i = 0; i < LOOP; ++i)
			AAT<map_by_unsorted_vector>(spA, spAAT);
		AAT<>(spA, spAAT);
	}
	convert(spAAT, AAT1);
	cout << norm(AAT0-AAT1) << endl;

	cout << "AA should be zero" << endl;
	csc<double> spB, spAB;
	trans(spA, spB);
#if ENABLE_DEPRECATED
	{timer t("spm_dmm");
		for(size_t i = 0; i < LOOP; ++i)
			spm_dmm(false, spA, false, spB, spAB);
	}
	convert(spAB, AAT1);
	cout << norm(AAT1-AAT0) << endl;
#endif
	{timer t("MM-sort");
		for(size_t i = 0; i < LOOP; ++i)
			MM<map_by_sorted_vector>(false, spA, false, spB, spAB);
	}
	{timer t("MM-map  ");
		for(size_t i = 0; i < LOOP; ++i)
			MM<std_map>(false, spA, false, spB, spAB);
	}
	{timer t("MM-unsort");
		for(size_t i = 0; i < LOOP; ++i)
			MM<map_by_unsorted_vector>(false, spA, false, spB, spAB);
	}
	convert(spAB, AAT1);
	cout << norm(AAT1-AAT0) << endl;

	cout << "add should be zero" << endl;
	matrix<double> B = floor(pow(rand<double>(m, k), 10)*100), C, D;
	csc<double> spC, spD;
	convert(B, spB, 1e-4);
	convert(spB, B);
	convert(spA, A);

	LOOP *= 20;
#if ENABLE_DEPRECATED
	{timer t("spm_mpm");
		for(size_t i = 0; i < LOOP; ++i)
			spm_mpm(spA, spB, spC);
	}
	convert(spC, C);
	cout << norm(A+B - C) << endl;
#endif
	csc_by_vm<double, ptrdiff_t, map_by_sorted_vector> vmE;
	{timer t("acc vm\t");
		for(size_t i = 0; i < LOOP; ++i) {
			vmE.resize(spA.size(1), spA.size(2));
			acc(spA, vmE);
			acc(spB, vmE);
			convert(vmE, spC);
		}
	}
	convert(spC, C);
	cout << norm(A+B - C) << endl;

	cout << (is_sorted_csc(spC)?"sorted":"unsorted") << endl;

	{timer t("acc csc\t");
		for(size_t i = 0; i < LOOP; ++i) {
			spC.val()(colon()) = 0;
			acc(spA, spC);
			acc(spB, spC);
		}
	}
	{timer t("acc csc sort");
		for(size_t i = 0; i < LOOP; ++i) {
			spC.val()(colon()) = 0;
			acc_sorted(spA, spC);
			acc_sorted(spB, spC);
		}
	}
	convert(spC, C);
	cout << norm(A+B - C) << endl;

	vmE.resize(spA.size(1), spA.size(2));
	acc(spA, vmE);
	acc(spB, vmE);
	convert(vmE, spC);
	convert(spC, C);
	cout << norm(A+B - C) << endl;

	cout << "check resize: nonzero then zero." << endl;
	cout << nnz(vmE) << endl;
	vmE.resize(spA.size(1), spA.size(2));
	cout << nnz(vmE) << endl;

	return 0;
}

int main(int argc, char *argv[])
{
//	if(test_basic_op())
//		return __LINE__;
//	test_solver();
//	test_boost();
	return 0;
}
