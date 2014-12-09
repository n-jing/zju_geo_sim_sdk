#ifndef PROBLEM_SIMPLE_FUNCTION_H
#define PROBLEM_SIMPLE_FUNCTION_H

#include <zjucad/matrix/io.h>
#include <jtflib/function/function.h>
#include <jtflib/function/func_aux.h>

#include "third_solver_ipopt.h"
#include "jtf_func_to_obj_func.h"

namespace zjucad
{
namespace LSCO
{

template <typename MATRIX_TYPE>
struct kernel_traits;

typedef jtf::function::functionN1_t<double, int32_t> jtf_func;
typedef std::shared_ptr<jtf_func> jtf_func_cons_ptr;
typedef std::vector<jtf_func_cons_ptr> func_container;

typedef hj::function::function_t<double, int32_t> hj_func;
typedef std::shared_ptr<const hj_func> hj_func_cons_ptr;

class one_func : public jtf_func
{
public :
    one_func(const size_t node_number,
             const size_t idx0,
             const size_t idx1)
        :node_number_(node_number), idx0_(idx0), idx1_(idx1) {}
    virtual ~one_func() {}
    virtual size_t dim() const
    {
        return node_number_;
    }
    virtual int val(const val_type *x, val_type &v) const
    {
        v += 0.5*(x[idx0_] - x[idx1_])*(x[idx0_]-x[idx1_]);
        return 0;
    }
    virtual int gra(const val_type *x, size_t &nnz, val_type *g, int_type *idx)
    {
        if(g == 0 && idx == 0) {
            nnz = 2;
            return 0;
        }

        g[0] = x[idx0_]-x[idx1_];
        idx[0] = idx0_;

        g[1] = x[idx1_]-x[idx0_];
        idx[1] = idx1_;

        return 0;
    }
    virtual int gra(const val_type *x, val_type *g)
    {
        g[idx0_] += x[idx0_]-x[idx1_];
        g[idx1_] += x[idx1_]-x[idx0_];
        return 0;
    }

    virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                    val_type *h, int_type *ptr, int_type *idx, val_type alpha)
    {
        if(h == 0  && idx == 0 && ptr == 0) {
            nnz = 4;
            format = 1;
            return 0;
        }
        if(h == 0 && idx != 0 && ptr != 0) {

            ptr[idx0_+1] = ptr[idx0_]+2;
            idx[ptr[idx0_]+0] = idx0_;
            idx[ptr[idx0_]+1] = idx1_;

            ptr[idx1_+1] = ptr[idx1_]+2;
            idx[ptr[idx1_]+0] = idx0_;
            idx[ptr[idx1_]+1] = idx1_;

            return 0;
        }
        if(h != 0 && idx != 0 && ptr != 0) {
            jtf::function::add_to_csc(h,ptr,idx, idx0_, idx0_,1.0);
            jtf::function::add_to_csc(h,ptr,idx, idx0_, idx1_,-1.0);

            jtf::function::add_to_csc(h,ptr,idx, idx1_, idx0_,-1.0);
            jtf::function::add_to_csc(h,ptr,idx, idx1_, idx1_,1.0);
        }
        return 0;
    }
    virtual int hes_block(const val_type *x, val_type *h, val_type alpha)
    {
        return __LINE__;
    }
private:
    const size_t node_number_;
    const size_t idx0_, idx1_;
};


class fix_func : public jtf_func
{
public :
    fix_func(const size_t node_number,
             const size_t idx0,
             const double v)
        :node_number_(node_number), idx0_(idx0), v_(v) {}
    virtual ~fix_func() {}
    virtual size_t dim() const
    {
        return node_number_;
    }
    virtual int val(const val_type *x, val_type &v) const
    {
        v += x[idx0_] - v_;
        return 0;
    }
    virtual int gra(const val_type *x, size_t &nnz, val_type *g, int_type *idx)
    {
        if(g == 0 && idx == 0) {
            nnz = 1;
            return 0;
        }

        g[0] = 1.0;
        idx[0] = idx0_;


        return 0;
    }
    virtual int gra(const val_type *x, val_type *g)
    {
        g[idx0_] += 1.0;

        return 0;
    }

    virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                    val_type *h, int_type *ptr, int_type *idx, val_type alpha)
    {
        nnz = 0;
        return 0;
    }
    virtual int hes_block(const val_type *x, val_type *h, val_type alpha)
    {
        return __LINE__;
    }
private:
    const size_t node_number_;
    const size_t idx0_;
    const double v_;
};


void test_simple()
{
    // y= 0.5*(x0-x1)^2
    // s.t. x0=0
    //      x1 > 2
    using namespace std;

    zjucad::matrix::matrix<double> test_node = zjucad::matrix::zeros<double>(2,1);
    jtf_func_cons_ptr total_obj(new one_func(test_node.size(),0,1));

    ////////////////////////////////////////////////////////////////////////////////////////////////

    typedef kernel_traits<typename hj::sparse::csc<double,int32_t> > kernel_type;

    std::shared_ptr<func_container > eqn_func(new func_container),
        ineqn_func(new func_container), test(new func_container);

    eqn_func->push_back(jtf_func_cons_ptr(new fix_func(test_node.size(), 0,1.0)));

    ineqn_func->push_back(jtf_func_cons_ptr(new fix_func(test_node.size(), 1,2.0)));

    third_solver_ipopt::func_ptr jftof(
        new jtf_func_to_obj_func<kernel_type>(total_obj, eqn_func, test));


    third_solver_ipopt tsb(jftof);
    boost::property_tree::ptree pt;
    pt.put("iter.value",10);
    tsb.solve(test_node, pt);
    std::cerr << test_node << std::endl;
}
}
}

#endif // PROBLEM_SIMPLE_FUNCTION_H
