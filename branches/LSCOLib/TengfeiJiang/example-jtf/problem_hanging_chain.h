#ifndef PROBLEM_HANGING_CHAIN_H
#define PROBLEM_HANGING_CHAIN_H

#include <jtflib/function/function.h>
#include <jtflib/function/func_aux.h>
#include <hjlib/sparse/sparse.h>
#include <zjucad/matrix/itr_matrix.h>

namespace zjucad
{
namespace LSCO
{
template <typename MATRIX_TYPE>
struct kernel_traits;

typedef jtf::function::functionN1_t<double,hj::sparse::idx_type> jtf_func;
typedef boost::shared_ptr<const jtf_func> jtf_func_cons_ptr;

class hanging_chain_func : public jtf_func
{
public:
    hanging_chain_func(const size_t segments,
                       const size_t idx)
        :segments_(segments), idx_(idx) {
        time_step_ = 1.0/segments;
    }
    virtual ~hanging_chain_func() {}

    virtual size_t dim() const {
        return 2 * (segments_-1);
    }

    virtual int val(const val_type *x, val_type &v) const {
        using namespace zjucad::matrix;
        itr_matrix<const val_type*> x0(segments_-1,2,x);
        v += time_step_ * x0(idx_,0) * std::sqrt(1+x0(idx_,1)*x0(idx_,1));

        return 0;
    }

    virtual int gra(const val_type *x, val_type *g) {
        using namespace zjucad::matrix;
        itr_matrix<const val_type*> x0(segments_-1,2,x);

        g[idx_] += time_step_ * std::sqrt(1+x0(idx_,1)*x0(idx_,1));
        g[idx_ + segments_-1] += time_step_ * x0(idx_,0) * x0(idx_,1)/(std::sqrt(1+x0(idx_,1)*x0(idx_,1)));

        return 0;
    }

    virtual int gra(const val_type *x, size_t &nnz, val_type *g, int_type *idx) {
        using namespace zjucad::matrix;
        itr_matrix<const val_type*> x0(segments_-1,2,x);

        if(g == 0 && idx == 0) {
            nnz = 2;
            return 0;
        }

        g[0] = time_step_ * std::sqrt(1+x0(idx_,1)*x0(idx_,1));
        g[1] = time_step_ * x0(idx_,0) * x0(idx_,1)/(std::sqrt(1+x0(idx_,1)*x0(idx_,1)));
        idx[0] = idx_;
        idx[1] = idx_+segments_-1;

        return 0;
    }

    virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                    val_type *h, int_type *ptr, int_type *idx, val_type alpha) {
        if(h == 0 && idx == 0 && ptr == 0) {
            format = 1;
            nnz = 3;
            return 0;
        }
        if(h == 0 && idx != 0 && ptr != 0) {
            ptr[idx_+1] = ptr[idx_] + 1;
            idx[ptr[idx_]] = idx_ + segments_-1;

            ptr[idx_+segments_-1] = ptr[idx_+1];
            ptr[idx_+segments_] = ptr[idx_+segments_-1] + 2;
            idx[ptr[idx_+segments_-1]+0] = idx_;
            idx[ptr[idx_+segments_-1]+1] = idx_+segments_-1;
        }
        if(h != 0 && idx != 0 && ptr != 0) {
            using namespace zjucad::matrix;
            itr_matrix<const val_type*> x0(segments_-1,2,x);

            const val_type sqrt_g = std::sqrt(1+x0(idx_,1)*x0(idx_,1));
            jtf::function::add_to_csc(h,ptr,idx, static_cast<int_type>(idx_+segments_-1), static_cast<int_type>(idx_), time_step_ * x0(idx_,1)/sqrt_g);
            jtf::function::add_to_csc(h,ptr,idx, static_cast<int_type>(idx_), static_cast<int_type>(idx_+segments_-1), time_step_ * x0(idx_,1)/sqrt_g);
            jtf::function::add_to_csc(h,ptr,idx, static_cast<int_type>(idx_+segments_-1), static_cast<int_type>(idx_+segments_-1), time_step_ * x0(idx_,0)/std::pow(sqrt_g,3));
        }
        return 0;
    }

    virtual int hes_block(const val_type *x, val_type *h, val_type alpha) {
        return __LINE__;
    }

private:
    const size_t segments_;
    double time_step_;
    const size_t idx_;
};

class length_constraints : public jtf_func
{
public:
    length_constraints(const size_t segments,
                       const double L):segments_(segments), L_(L) {
        time_step_ = 1.0/segments;
    }
    virtual ~length_constraints() {}

    virtual size_t dim() const {
        return 2*(segments_-1);
    }
    virtual int val(const val_type *x, val_type &v) const {
        using namespace zjucad::matrix;
        itr_matrix<const val_type*> x0(segments_-1,2,x);
        v -= L_;
        for(size_t i = 0; i < segments_-1; ++i) {
            v += time_step_ *std::sqrt(1+x0(i,1)*x0(i,1));
        }
        return 0;
    }
    virtual int gra(const val_type *x, val_type *g) {
        using namespace zjucad::matrix;
        itr_matrix<const val_type*> x0(segments_-1,2,x);
        for(size_t i = 0; i < segments_-1; ++i) {
            g[i+segments_-1] += time_step_ * x0(i,1)/std::sqrt(1+x0(i,1)*x0(i,1));
        }
        return 0;
    }
    virtual int gra(const val_type *x, size_t &nnz, val_type *g, int_type *idx) {
        if(g == 0 && idx == 0) {
            nnz = segments_-1;
            return 0;
        }

        using namespace zjucad::matrix;
        itr_matrix<const val_type*> x0(segments_-1,2,x);
        for(size_t i = 0; i < segments_-1; ++i) {
            g[i] = time_step_*x0(i,1)/std::sqrt(1+x0(i,1)*x0(i,1));
            idx[i] = i+segments_-1;
        }

        return 0;
    }
    virtual int hes(const val_type *x, size_t &nnz, size_t &format, val_type *h,
                    int_type *ptr, int_type *idx, val_type alpha) {
        if(h == 0 && idx == 0 && ptr == 0) {
            nnz = segments_-1;
            format = 1;
            return 0;
        }
        if(h == 0 && idx != 0 && ptr != 0) {
            for(size_t i = 0 ; i < segments_-1; ++i) {
                ptr[i+segments_+1] = ptr[i+segments_] + 1;
                idx[ptr[i+segments_]] = i + segments_;
            }
        }
        if(h != 0 && idx != 0 && ptr != 0) {
            using namespace zjucad::matrix;
            itr_matrix<const val_type*> x0(segments_-1,2,x);
            for(size_t i  = 0 ; i < segments_ -1; ++i) {
                jtf::function::add_to_csc(h,ptr,idx, static_cast<int_type>(i+segments_-1), static_cast<int_type>(i+segments_-1), time_step_ * std::pow(1+x0(i,1)*x0(i,1),-1.5));
            }
        }
        return 0;
    }
    virtual int hes_block(const val_type *x, val_type *h, val_type alpha) {
        return __LINE__;
    }
private:
    const size_t segments_;
    const double L_;
    double time_step_;
};

class derivative_constraint : public jtf_func
{
    // g_i = (x_i+1 - xi-1)/2*time_step
public:
    derivative_constraint(const size_t segments, const size_t idx,
                          const double a, const double b)
        : segments_(segments), time_step_(1.0/segments), idx_(idx), a_(a), b_(b) {}
    virtual ~derivative_constraint() {}
    virtual size_t dim() const {
        return 2 * (segments_-1);
    }
    virtual int val(const val_type *x, val_type &v) const {
        using namespace zjucad::matrix;
        itr_matrix<const val_type*> x0(segments_-1,2,x);
        if(idx_ == 0) {
            v += x0(idx_,1) - (x0(idx_+1,0) -a_)/(2*time_step_);
        } else if(idx_ == segments_-2) {
            v += x0(idx_,1) - (b_-x0(idx_-1,0))/(2*time_step_);
        } else {
            v += x0(idx_,1) - (x0(idx_+1,0)-x0(idx_-1,0))/(2*time_step_);
        }

        return 0;
    }

    virtual int gra(const val_type *x, val_type *g) {
        if(idx_ == 0) {
            g[idx_+segments_-1] += 1.0;
            g[idx_+1] += -1.0/(2*time_step_);
        } else if(idx_ == segments_-2) {
            g[idx_+segments_-1] += 1.0;
            g[idx_-1] += 1.0/(2*time_step_);
        } else {
            g[idx_+segments_-1] += 1.0;
            g[idx_-1] += 1.0/(2*time_step_);
            g[idx_+1] += -1.0/(2*time_step_);
        }
        return 0;
    }

    virtual int gra(const val_type *x, size_t &nnz, val_type *g, int_type *idx) {
        if(g == 0 && idx == 0) {
            if(idx_ == 0 || idx_ == segments_-2) nnz = 2;
            else
                nnz = 3;
            return 0;
        }

        if(idx_ == 0) {
            g[0] = 1.0;
            idx[0] = idx_+segments_-1;
            g[1] = -1.0/(2*time_step_);
            idx[1] = idx_+1;
        } else if(idx_ == segments_-2) {
            g[0] = 1.0;
            idx[0] = idx_+segments_-1;
            g[1] = 1.0/(2*time_step_);
            idx[1] = idx_-1;
        } else {
            g[0] = 1.0;
            idx[0] = idx_+segments_-1;
            g[1] = 1.0/(2*time_step_);
            idx[1] = idx_-1;
            g[2] = -1.0/(2*time_step_);
            idx[2] = idx_+1;
        }
        return 0;
    }

    virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                    val_type *h, int_type *ptr, int_type *idx, val_type alpha) {
        nnz = 0;
        return 0;
    }

    virtual int hes_block(const val_type *x, val_type *h, val_type alpha) {
        return __LINE__;
    }
private:
    const size_t segments_;
    const double time_step_;
    const size_t idx_;
    const double a_, b_;
};

void hanging_chain(const double a, const double b, const double L, const size_t segments)
{
    typedef std::vector<jtf_func_cons_ptr> jtf_func_container;

    zjucad::matrix::matrix<double> xg = zjucad::matrix::zeros<double>(segments-1,2);

    jtf_func_container jfc;
    {
        for(size_t i = 0; i < segments-1; ++i) {
            jfc.push_back(jtf_func_cons_ptr(new hanging_chain_func(segments,i)));
        }
    }

    jtf_func_cons_ptr obj(new jtf::function::sum_function<double,jtf_func::int_type,jtf::function::SMART_BOOST>(jfc));

    boost::shared_ptr<jtf_func_container> eqn_constraints(new jtf_func_container);
    {
        eqn_constraints->push_back(jtf_func_cons_ptr(new length_constraints(segments,L)));
        for(size_t i = 0; i < segments-1; ++i) {
            eqn_constraints->push_back(jtf_func_cons_ptr(new derivative_constraint(segments,i,a,b)));
        }
    }

    typedef kernel_traits<typename hj::sparse::csc<double, hj::sparse::idx_type> > kernel_type;

    std::cerr << jtf::function::gra_err(*obj, &xg[0]) << std::endl;
    std::cerr << jtf::function::hes_err(*obj, &xg[0]) << std::endl;
}
}
}
#endif // PROBLEM_HANGING_CHAIN_H
