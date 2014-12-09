#ifndef MESH_FUNCTION_H
#define MESH_FUNCTION_H

#include <hjlib/function/function.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>

///////////////////////////////////////////////////
/////////// low level func ////////////////////////

class fix_node : public hj::function::function_t<double,int32_t>
{
public:
    fix_node(const zjucad::matrix::matrix<double> & position,
             const size_t & node_num,
             const size_t & p,
             const double w)
        : position_(position), node_num_(node_num), p_(p), w_(w){}

    virtual size_t dim_of_x(void) const {
        return node_num_*3;
    }
    virtual size_t dim_of_f(void) const {
        return 3;
    }
    virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
        const zjucad::matrix::itr_matrix<const double*> T(3, node_num_, x);
        zjucad::matrix::itr_matrix<double*> f0(3,1,f);
        f0  = T(zjucad::matrix::colon(),p_) - position_;
        f0 *= w_;
        return 0;
    }
    virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                    int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
        for(int fi = 0; fi < 3; ++fi) {
            ptr[fi+1] = ptr[fi]+1;
            idx[ptr[fi]] = 3 * p_ + fi;
            val[ptr[fi]] = w_;
        }
        return 0;
    }
    virtual size_t jac_nnz(void) const {
        return 1*dim_of_f();
    }
private:
    const zjucad::matrix::matrix<double> position_;
    const size_t node_num_;
    const size_t p_;
    const double w_;
};


class align_orthog_dir_func : public hj::function::function_t<double, int32_t>
{
public:
    align_orthog_dir_func(const zjucad::matrix::matrix<double> & node,
                          const size_t & node_idx,
                          const zjucad::matrix::matrix<double> & orth_dir,
                          const double w)
        : position_(node(zjucad::matrix::colon(),node_idx)),
          node_num_(node.size(2)), orthog_dir_(orth_dir),  p_(node_idx), w_(w){}

    virtual size_t dim_of_x(void) const {
        return node_num_*3;
    }
    virtual size_t dim_of_f(void) const {
        return 1;
    }
    virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
        const zjucad::matrix::itr_matrix<const double*> T(3, node_num_, x);
        zjucad::matrix::matrix<double> displace = T(zjucad::matrix::colon(), p_) - position_;
        *f = dot(orthog_dir_, displace);
        *f *= w_;
        return 0;
    }
    virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                    int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
        ptr[1] = ptr[0]+3;
        for(int fi = 0; fi < 3; ++fi) {

            idx[ptr[0] + fi] = 3 * p_ + fi;
            val[ptr[0] + fi] = w_ * orthog_dir_[fi];
        }
        return 0;
    }
    virtual size_t jac_nnz(void) const {
        return 3;
    }
private:
    const zjucad::matrix::matrix<double> position_;
    const zjucad::matrix::matrix<double> orthog_dir_;
    const size_t node_num_;
    const size_t p_;
    const double w_;
};

#endif //MESH_FUNCTION_H
