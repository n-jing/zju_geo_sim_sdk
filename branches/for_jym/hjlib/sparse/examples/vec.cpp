#include "sparse.h"
#include "io.h"

#include <boost/numeric/ublas/vector_sparse.hpp>

#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>

using namespace std;
using namespace boost::numeric;
using namespace hj::sparse;
template <typename T>
void test_sparse_vec(int LEN, const vector<ptrdiff_t> &idx, int LOOP, int nnz)
{
	for(int itr = 0; itr < LOOP; ++itr) {
		T v(LEN, nnz);
		for(size_t i = 0; i < idx.size(); ++i)
			v[idx[i]] += 1;
	}
}

void test_performance(const vector<ptrdiff_t> &idx)
{
	clock_t beg, elapse;

	const int NNZ = 100; //*max_element(idx.begin(), idx.end());
	const int LEN = 10000, LOOP = 100000000L/NNZ/idx.size();

	double time_unit;
	for(int i = 0; i < 1; ++i) {
		beg = clock();
		test_sparse_vec<map_vec<> >(LEN, idx, LOOP, NNZ);
		elapse = clock()-beg;
		time_unit = elapse;
		cout << "map: " << elapse/time_unit << endl;

		beg = clock();
		test_sparse_vec<ublas::mapped_vector<double, map<ptrdiff_t, double> > >(LEN, idx, LOOP, NNZ);
		elapse = clock()-beg;
		cout << "ublas: " << elapse/time_unit << endl;


		beg = clock();
		test_sparse_vec<map_vec<double, ptrdiff_t, map_by_unsorted_vector> >(LEN, idx, LOOP, NNZ);
		elapse = clock()-beg;
		cout << "unsort: " << elapse/time_unit << endl;

		beg = clock();
		test_sparse_vec<map_vec<double, ptrdiff_t, map_by_sorted_vector> >(LEN, idx, LOOP, NNZ);
		elapse = clock()-beg;
		cout << "sort: " << elapse/time_unit << endl;

		cout << endl;
	}
}

void test_insert(void)
{
	const int NNZ = 50;
	vector<ptrdiff_t> idx(NNZ);

	for(int i = 0; i < 5; ++i) {
		int NUM = 2 << i;
		idx.resize(NNZ*(NUM+1));
		for(size_t i = 0; i < idx.size(); ++i)
			idx[i] = i;
		random_shuffle(idx.begin(), idx.end());

		test_performance(idx);
	}
}


template <typename T>
void check_correct(void)
{
	T v(10);
//	cout << v;  // STRANGE: this code can prevent the warning: coo.h:86: warning: dereferencing pointer ‘<anonymous>’ does break strict-aliasing rules
	v[0] = 20;
	v[1] = 2;
	v[2] = -0.2;
	v[3] = 0.3;
	v[4] += 6;
	v[2] += 0.7;
	v[2] += 4;
	cout << nnz(v) << endl << v << endl;
}

void check_correct()
{
	check_correct<map_vec<double, ptrdiff_t, map_by_unsorted_vector> >();
	check_correct<map_vec<double, ptrdiff_t, map_by_sorted_vector> >();
	check_correct<map_vec<double, ptrdiff_t> >();
}

int main(void)
{
//	test_insert();
//	test_conflict();
	check_correct();
	return 0;
}
