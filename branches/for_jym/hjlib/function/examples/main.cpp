#include <stdint.h>
#include <cmath>
#include <vector>
#include <iterator>
#include <iostream>
#include <memory>
#include <cassert>

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>

#include <boost/shared_ptr.hpp>

#include "../include/function.h"
#include "../include/func_aux.h"

using namespace std;
using namespace zjucad::matrix;
using namespace hj::function;
using namespace boost;

#define USE_VALUE_TYPE double
#define USE_INT_TYPE int32_t

class fun0 : public function_t<USE_VALUE_TYPE, USE_INT_TYPE>
{
public:

	static size_t s_dim_of_x(void) { return 4; }
	static size_t s_dim_of_f(void) { return 2; }

	virtual size_t dim_of_x(void) const { return s_dim_of_x(); }
	virtual size_t dim_of_f(void) const { return s_dim_of_f(); }

	virtual int val(const value_type *x, value_type *f, func_ctx *ctx = 0) const {
		f[0] = cos(x[0]);
		f[1] = x[1]*x[2];

		return 0;
	}
	virtual int jac(const value_type *x, value_type *val, int_type *ptr = 0, int_type *idx = 0, func_ctx *ctx = 0) const {
		ptr[1] = ptr[0]+1;
		val[ptr[0]] = -sin(x[0]);
		idx[ptr[0]] = 0;

		ptr[2] = ptr[1]+2;
		val[ptr[1]] = x[2];
		idx[ptr[1]] = 1;
		val[ptr[1]+1] = x[1];
		idx[ptr[1]+1] = 2;

		return 0;
	}

	virtual size_t jac_nnz(void) const {
		return 1+2;
	}
};

class fun1 : public function_t<USE_VALUE_TYPE, USE_INT_TYPE>
{
public:
	virtual size_t dim_of_x(void) const { return 4; }
	virtual size_t dim_of_f(void) const { return 3; }

	virtual int val(const value_type *x, value_type *f, func_ctx *ctx = 0) const {
		f[0] = 2-x[1];
		f[1] = x[3]-x[0];
		f[2] = x[2]-sin(x[0]);

		return 0;
	}
	virtual int jac(const value_type *x, value_type *val, int_type *ptr = 0, int_type *idx = 0, func_ctx *ctx = 0) const {
		ptr[1] = ptr[0]+1;
		val[ptr[0]] = -1;
		idx[ptr[0]] = 1;

		ptr[2] = ptr[1]+2;
		val[ptr[1]] = 1;
		idx[ptr[1]] = 3;
		val[ptr[1]+1] = -1;
		idx[ptr[1]+1] = 0;

		ptr[3] = ptr[2]+2;
		val[ptr[2]] = 1;
		idx[ptr[2]] = 2;
		val[ptr[2]+1] = -cos(x[0]);
		idx[ptr[2]+1] = 0;

		return 0;
	}

	virtual size_t jac_nnz(void) const {
		return 1+2+2;
	}
};

template <typename OS, typename T>
OS &operator << (OS &os, const vector<T> &v)
{
	copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
	return os;
}

void usage0(const vector<fun0::value_type> &x)
{
	fun0 f0;
	fun1 f1;

	assert(f0.dim_of_x() == f1.dim_of_x());

	// checking val
	cout << "\nbegin to check val." << endl;

	vector<fun0::value_type> v0(f0.dim_of_f());
	f0.val(&x[0], &v0[0]);
	cout << v0 << endl;
	
	auto_ptr<hj::function::function> sf0(f0/10.0);
	vector<fun0::value_type> sv0(sf0->dim_of_f());
	sf0->val(&x[0], &sv0[0]);
	cout << sv0 << endl;

	vector<fun1::value_type> v1(f1.dim_of_f());
	f1.val(&x[0], &v1[0]);
	cout << v1 << endl;

	vector<hj::function::function *> cf0(2);
	cf0[0] = sf0.get();
	cf0[1] = &f1;
	auto_ptr<hj::function::function> cf(new_catenated_function<fun1::value_type, fun1::int_type>(cf0));

	vector<fun0::value_type> cfv(cf->dim_of_f());
	cf->val(&x[0], &cfv[0]);
	cout << "check f: " << cfv << endl;

	// checking jac
	cout << "\nbegin to check jac." << endl;

	vector<fun0::value_type> val(cf->jac_nnz());
	vector<fun0::int_type> ptr(cf->dim_of_f()+1), idx(cf->jac_nnz());

	cf->jac(&x[0], &val[0], &ptr[0], &idx[0]);
	cout << ptr << "\n" << val << "\n" << idx << endl;

	matrix<fun0::value_type> jac0(cf->dim_of_f(), cf->dim_of_x()),
		jac1(cf->dim_of_f(), cf->dim_of_x());

	for(size_t ci = 0; ci < ptr.size()-1; ++ci) {
		for(size_t nzi = ptr[ci]; nzi < ptr[ci+1]; ++nzi) {
			jac0(ci, idx[nzi]) = val[nzi];
		}
	}
	calc_jac(*cf, &x[0], &jac1[0]);
	cout << "jac max err: " << max(fabs(jac1-jac0)) << endl;
}

hj::function::function *new_func(void)
{
	typedef boost::shared_ptr<const function_t<fun0::value_type, fun0::int_type> > const_func_ptr;
	boost::shared_ptr<vector<const_func_ptr> > cfs(new vector<const_func_ptr>);
	cfs->push_back(const_func_ptr(
					   const_func_ptr(new fun0)/10.0));
	cfs->push_back(const_func_ptr(new fun1));
	return new_catenated_function<fun0::value_type, fun0::int_type>(cfs);
}

void usage1(const vector<fun0::value_type> &x)
{
	boost::shared_ptr<hj::function::function> cf(new_func());
	size_t len = cf->dim_of_f();
	vector<fun0::value_type> cfv(len);
	cf->val(&x[0], &cfv[0]);
	cout << "check f: " << cfv << endl;
}

void type_test(void)
{
	cout << "\nbegin to check type." << endl;

	fun0 f0;

	cout << dynamic_cast<function_t<double, int32_t> *>(&f0) << endl;
	cout << dynamic_cast<function_t<float, int32_t> *>(&f0) << endl;
	cout << dynamic_cast<function_t<double, int64_t> *>(&f0) << endl;

	cout << f0.get_value_type() << " " << f0.get_int_type() << endl;

	const hj::function::function *f = &f0;
	cout << "may assert or crash due to improper type in debug." << endl;
	vector<float> x(f0.dim_of_x()), val(f0.dim_of_f());
	f->val(&x[0], &val[0]);
}

int main(int argc, char *argv[])
{
	vector<fun0::value_type> x(fun0::s_dim_of_x());
	for(size_t i = 0; i < x.size(); ++i)
		x[i] = rand()/static_cast<double>(RAND_MAX);

	cout << x << endl;

	usage0(x);
	usage1(x);

//	type_test();

	return 0;
}
