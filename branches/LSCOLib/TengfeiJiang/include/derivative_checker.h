#ifndef DERIVATIVE_CHECKER_H
#define DERIVATIVE_CHECKER_H

#include "kernel.h"
#include "objective_function.h"
#include <boost/static_assert.hpp>
#include <hjlib/sparse/sparse.h>
namespace zjucad
{
namespace LSCO
{
template <typename MATRIX_TYPE>
struct kernel_traits;

template <typename KERNEL_TYPE>
class derivative_checker
{
public:
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;

    BOOST_STATIC_ASSERT((boost::is_same<sparse_matrix_type, hj::sparse::csc<value_type> >::value));

    explicit derivative_checker(const objective_function<KERNEL_TYPE> &of):m_of(of) {}
    explicit derivative_checker(boost::shared_ptr<objective_function<KERNEL_TYPE> > &ptr):m_of(*ptr) {}

    template <typename OS>
    void check(OS &os, vector_type &x, const size_t order) const {
        switch(order) {
        case 1 : {
            check_first_order(os, x);
            break;
        }
        case 2 : {
            check_second_order(os, x);
            break;
        }
        case 3 : {
            check_first_order(os, x);
            check_second_order(os, x);
            break;
        }
        default:
            std::cerr << "# [error] invalid order, [1:first/2:second/3:all]" << std::endl;
            break;
        }
    }
protected:
    template <typename OS>
    void check_first_order(OS &os, vector_type &x)const {
        const size_t dim = m_of.nr_inputs();
        value_type val = 0;
        vector_type gra;
        KERNEL_TYPE::resize(dim, gra);
        KERNEL_TYPE::zero(gra);
        val = m_of.eval_func(x,gra);

        const value_type eps = 1e-6;
        value_type v[2]  = {0,0};
        for(size_t xi = 0; xi < dim; ++xi) {
            const value_type save = KERNEL_TYPE::get(xi,x);
            KERNEL_TYPE::set(xi, save-eps ,x);
            v[0] = m_of.eval_func(x, gra, false);
            KERNEL_TYPE::set(xi, save+eps ,x);
            v[1] = m_of.eval_func(x, gra, false);
            KERNEL_TYPE::set(xi, KERNEL_TYPE::get(xi,gra) - (v[1]-v[0])/(2*eps), gra);
            KERNEL_TYPE::set(xi, save, x);
        }
        os << "gra diff norm " << KERNEL_TYPE::nrm2(gra) << std::endl;
        os << "gra amax diff " << KERNEL_TYPE::amax(gra) << std::endl;
    }

    template <typename OS>
    void check_second_order(OS &os, vector_type &x)const {
        const size_t dim = m_of.nr_inputs();
        sparse_matrix_type H;
        m_of.eval_func_hessian(x,H,false);
        vector_type ga,gb;
        zjucad::matrix::matrix<value_type> H_m(dim,dim);
        hj::sparse::convert(H,H_m);

        KERNEL_TYPE::resize(x.size(), ga);
        KERNEL_TYPE::resize(x.size(), gb);

        const value_type eps = 1e-6;
        for(size_t xi = 0; xi < dim; ++xi) {
            const value_type save = KERNEL_TYPE::get(xi,x);
            KERNEL_TYPE::set(xi, save+eps, x);
            KERNEL_TYPE::zero(ga);
            m_of.eval_func(x,ga, true);
            KERNEL_TYPE::set(xi,save-eps, x);
            KERNEL_TYPE::zero(gb);
            m_of.eval_func(x, gb, true);
            KERNEL_TYPE::set(xi,save, x);
            H_m(zjucad::matrix::colon(),xi) -= (ga-gb)/(2*eps);
        }
        os << "hes diff norm " << zjucad::matrix::norm(H_m) << std::endl;
        os << "hes amax diff " << zjucad::matrix::max(zjucad::matrix::fabs(H_m)) << std::endl;
    }

private:
    const objective_function<KERNEL_TYPE> &m_of;
    derivative_checker() {}
};
}
}
#endif // DERIVATIVE_CHECKER_H
