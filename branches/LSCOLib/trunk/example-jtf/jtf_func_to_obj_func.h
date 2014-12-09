#ifndef JTF_FUNCTION_TO_OBJ_FUNC_H
#define JTF_FUNCTION_TO_OBJ_FUNC_H

#include <jtflib/function/function.h>
#include "objective_function.h"
#include "kernel_zjucad_matrix.h"

namespace zjucad
{
namespace LSCO
{
template <typename MATRIX_TYPE>
struct kernel_traits;

template <typename KERNEL_TYPE>
class jtf_func_to_obj_func : public objective_function<KERNEL_TYPE>
{
public:
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;

    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef typename KERNEL_TYPE::int_type int_type;

    typedef jtf::function::functionN1_t<value_type, int_type> JTF_FUNC;
    typedef std::vector<std::shared_ptr<JTF_FUNC> > func_container;

    jtf_func_to_obj_func(std::shared_ptr<JTF_FUNC> jtf_obj_func,
                         std::shared_ptr<func_container> eqn_container,
                         std::shared_ptr<func_container> ineqn_constainer)
        :jtf_func_(jtf_obj_func), nnz_h(-1), format_h(-1) ,eqn_container_(eqn_container),
         ineqn_container_(ineqn_constainer)
    {
        if(!eqn_container_.get()) eqn_container_.reset(new func_container);
        if(!ineqn_container_.get()) ineqn_container_.reset(new func_container);
        static vector_type fake_x;
        KERNEL_TYPE::resize(nr_inputs(),fake_x);
        assemble(fake_x);
    }
    virtual ~jtf_func_to_obj_func() {}

    ///
    /// \brief nr_inputs return dimension of x
    /// \return
    ///
    virtual size_t nr_inputs()const
    {
        return jtf_func_->dim();
    }

    virtual size_t nr_equals()const
    {
        if(!eqn_container_.get()) return 0;
        return eqn_container_->size();
    }

    virtual size_t nr_inequals()const
    {
        if(!ineqn_container_.get()) return 0;
        return ineqn_container_->size();
    }

    virtual value_type eval(const vector_type& x,bool same_x)
    {
        // ignore same_x
        value_type f = 0;
        jtf_func_->val(&x[0], f);
        return f;
    }

    virtual void eval_func_gradient(const vector_type &x, vector_type &gradient, bool same_x)
    {
        KERNEL_TYPE::zero(gradient);
        jtf_func_->gra(&x[0],&gradient[0]);
    }

    virtual void eval_lagrangian_hessian(const vector_type& x,
                                         const vector_type& lambda,
                                         value_type coef,
                                         typename KERNEL_TYPE::matrix_builder& hessian,
                                         bool same_x)
    {
        // ignore same_x
        size_t format = 1;
        H_.val() *= 0;
        jtf_func_->hes(&x[0], nnz_h, format_h, &H_.val()[0], &H_.ptr()[0], &H_.idx()[0]);
        H_.val() *= coef;

        for(size_t eqi = 0; eqi < eqn_container_->size(); ++eqi) {
            if(eqn_H_[eqi].size(1) != 0 && eqn_H_[eqi].size(2) != -1) {
                eqn_H_[eqi].val() *= 0;
                size_t nnz = hj::sparse::nnz(eqn_H_[eqi]);
                if(nnz > 0)
                    (*eqn_container_)[eqi]->hes(&x[0], nnz , format, 0,
                                                &eqn_H_[eqi].ptr()[0], &eqn_H_[eqi].idx()[0]);
            }
        }

        for(size_t eqi = 0; eqi < ineqn_container_->size(); ++eqi) {
            if(ineqn_H_[eqi].size(1) != 0 && ineqn_H_[eqi].size(2) != -1) {
                ineqn_H_[eqi].val() *= 0;
                size_t nnz = hj::sparse::nnz(ineqn_H_[eqi]);
                if(nnz > 0)
                    (*ineqn_container_)[eqi]->hes(&x[0], nnz , format, &ineqn_H_[eqi].val()[0],
                                                  &ineqn_H_[eqi].ptr()[0], &ineqn_H_[eqi].idx()[0]);
            }
        }

        KERNEL_TYPE::clear(hessian);
        if(!ij2double_ptr_hes_.empty()) {
            //for(const auto & term : ij2double_ptr_hes_) {
            for(typename std::map<std::pair<int_type,int_type>, std::vector<std::pair<int_type,value_type*> > >::const_iterator
                    cit = ij2double_ptr_hes_.begin(); cit != ij2double_ptr_hes_.end(); ++cit) {
                const std::pair<std::pair<int_type,int_type>, std::vector<std::pair<int_type,value_type*> > > & term = *cit;
                value_type h_temp = 0;
                // for(const auto & ptr : term.second) {
                for(typename std::vector<std::pair<int_type,value_type*> >::const_iterator pit = term.second.begin();
                        pit != term.second.end(); ++pit) {
                    const std::pair<int_type,value_type*>  & ptr = *pit;
                    double v = *ptr.second;
                    if(ptr.first != -1)// if is constraint
                        v *= lambda[ptr.first];
                    h_temp += v;
                }

                hessian.push_back(std::make_pair(std::make_pair(term.first.first, term.first.second), h_temp));
            }
        }
    }
    virtual void eval_constraint_value(const vector_type &x, vector_type &c, bool same_x)
    {
        // ignore same_x
        size_t i = 0;
        value_type f = 0;
        for(; i < eqn_container_->size(); ++i) {
            f = 0;
            (*eqn_container_)[i]->val(&x[0], f);
            c[i] = f;
        }
        for(; i < eqn_container_->size() + ineqn_container_->size(); ++i) {
            f = 0;
            (*ineqn_container_)[i-eqn_container_->size()]->val(&x[0], f);
            c[i] = f;
        }
    }
    virtual void eval_constraint_gradient(const vector_type& x,
                                          typename KERNEL_TYPE::matrix_builder& c_gradient,
                                          bool same_x)
    {
        // ignore same x
        static std::vector<value_type> g;
        static std::vector<int_type> idx;
        KERNEL_TYPE::clear(c_gradient);
        size_t ci = 0, nnz;
        for(; ci < eqn_container_->size(); ++ci) {
            (*eqn_container_)[ci]->gra(&x[0], nnz, 0,0);
            g.resize(nnz);
            idx.resize(nnz);
            (*eqn_container_)[ci]->gra(&x[0], nnz, &g[0], &idx[0]);
            for(size_t i = 0; i < nnz; ++i) {
                c_gradient.push_back(std::make_pair(std::make_pair(ci, idx[i]),g[i]));
            }
        }
        for(; ci < eqn_container_->size() + ineqn_container_->size(); ++ci) {
            (*ineqn_container_)[ci-eqn_container_->size()]->gra(&x[0], nnz, 0,0);
            g.resize(nnz);
            idx.resize(nnz);
            (*ineqn_container_)[ci-eqn_container_->size()]->gra(&x[0], nnz, &g[0], &idx[0]);
            for(size_t i = 0; i < nnz; ++i) {
                c_gradient.push_back(std::make_pair(std::make_pair(ci, idx[i]),g[i]));
            }
        }
    }

private:
    void assemble(const vector_type &x)
    {
        {
            if(H_.size(1) == 0 || H_.size(2) == 0) {
                jtf_func_->hes(&x[0], nnz_h, format_h,0,0,0);
                if(nnz_h > 0) {
                    H_.resize(jtf_func_->dim(), jtf_func_->dim(), nnz_h);
                    jtf_func_->hes(&x[0], nnz_h, format_h, 0, &H_.ptr()[0], &H_.idx()[0]);
                }
            }

            if(!eqn_container_->empty()) {
                eqn_H_.resize(eqn_container_->size());
                size_t nnz , format;
                for(size_t eqi = 0; eqi < eqn_container_->size(); ++eqi) {
                    (*eqn_container_)[eqi]->hes(&x[0], nnz, format, 0,0,0);
                    if(nnz > 0) {
                        eqn_H_[eqi].resize(nr_inputs(),nr_inputs(), nnz);
                        (*eqn_container_)[eqi]->hes(&x[0],nnz,format, 0, &eqn_H_[eqi].ptr()[0],
                                                    &eqn_H_[eqi].idx()[0]);
                    }
                }
            }

            if(!ineqn_container_->empty()) {
                ineqn_H_.resize(ineqn_container_->size());
                size_t nnz , format;
                for(size_t eqi = 0; eqi < ineqn_container_->size(); ++eqi) {
                    (*ineqn_container_)[eqi]->hes(&x[0], nnz, format, 0,0,0);
                    if(nnz > 0) {
                        ineqn_H_[eqi].resize(nr_inputs(),nr_inputs(), nnz);
                        (*ineqn_container_)[eqi]->hes(&x[0],nnz,format, 0, &ineqn_H_[eqi].ptr()[0],
                                                      &ineqn_H_[eqi].idx()[0]);
                    }
                }
            }

            if(ij2double_ptr_hes_.empty()) {
                if(H_.size(1) != 0 || H_.size(2) != 0) {
                    for(size_t i = 0; i < H_.ptr().size() - 1; ++i) {
                        for(size_t pi = H_.ptr()[i]; pi < H_.ptr()[i+1]; ++pi) {
                            ij2double_ptr_hes_[std::make_pair(i,H_.idx()[pi])].push_back(
                                std::make_pair(-1,&H_.val()[pi]));
                        }
                    }
                }
                for(size_t eqi = 0; eqi < eqn_H_.size(); ++eqi) {
                    if(eqn_H_[eqi].size(1) == 0 || eqn_H_[eqi].size(2) == 0) continue;
                    for(size_t i = 0; i < eqn_H_[eqi].ptr().size() - 1; ++i) {
                        for(size_t pi = eqn_H_[eqi].ptr()[i]; pi < eqn_H_[eqi].ptr()[i+1]; ++pi) {
                            ij2double_ptr_hes_[std::make_pair(i,eqn_H_[eqi].idx()[pi])].push_back(
                                std::make_pair(eqi,&eqn_H_[eqi].val()[pi]));
                        }
                    }
                }
                for(size_t eqi = 0; eqi < ineqn_H_.size(); ++eqi) {
                    if(ineqn_H_[eqi].size(1) == 0 || ineqn_H_[eqi].size(2) == 0) continue;
                    for(size_t i = 0; i < ineqn_H_[eqi].ptr().size() - 1; ++i) {
                        for(size_t pi = ineqn_H_[eqi].ptr()[i]; pi < ineqn_H_[eqi].ptr()[i+1]; ++pi) {
                            ij2double_ptr_hes_[std::make_pair(i,ineqn_H_[eqi].idx()[pi])].push_back(
                                std::make_pair(eqi+eqn_container_->size(),&ineqn_H_[eqi].val()[pi]));
                        }
                    }
                }
            }
        }
    }
private:
    std::shared_ptr<JTF_FUNC> jtf_func_;
    std::shared_ptr<func_container> eqn_container_, ineqn_container_;
    hj::sparse::csc<value_type,int_type> H_;
    std::vector<hj::sparse::csc<value_type, int_type> > eqn_H_, ineqn_H_;
    std::map<std::pair<int_type,int_type>, std::vector<value_type*> >  ij2double_ptr_g_jac_;
    std::map<std::pair<int_type,int_type>, std::vector<std::pair<int_type,value_type*> > > ij2double_ptr_hes_;
    size_t nnz_h, format_h;
};

template <typename KERNEL_TYPE>
class jtf_func_to_obj_func_L_D_BFGS : public objective_function_L_D_BFGS<KERNEL_TYPE>
{
public:
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;

    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef typename KERNEL_TYPE::int_type int_type;

    typedef jtf::function::functionN1_t<value_type, int_type> JTF_FUNC;
    typedef std::vector<std::shared_ptr<const JTF_FUNC> > func_container;

    jtf_func_to_obj_func_L_D_BFGS(std::shared_ptr<JTF_FUNC> jtf_obj_func,
                                  std::shared_ptr<func_container> eqn_container,
                                  std::shared_ptr<func_container> ineqn_constainer)
        :jtf_func_(jtf_obj_func), nnz_h(-1), format_h(-1) ,eqn_container_(eqn_container),
         ineqn_container_(ineqn_constainer)
    {
        if(!eqn_container_.get()) eqn_container_.reset(new func_container);
        if(!ineqn_container_.get()) ineqn_container_.reset(new func_container);
        static vector_type fake_x;
        KERNEL_TYPE::resize(nr_inputs(),fake_x);
        assemble(fake_x);
    }
    virtual ~jtf_func_to_obj_func_L_D_BFGS() {}

    virtual size_t nr_inputs()const
    {
        return jtf_func_->dim();
    }

    virtual size_t nr_equals()const
    {
        if(!eqn_container_.get()) return 0;
        return eqn_container_->size();
    }

    virtual size_t nr_inequals()const
    {
        if(!ineqn_container_.get()) return 0;
        return ineqn_container_->size();
    }

    virtual value_type eval(const vector_type& x,bool same_x)
    {
        // ignore same_x
        value_type f = 0;
        jtf_func_->val(&x[0], f);
        return f;
    }

    virtual void eval_func_gradient(const vector_type &x, vector_type &gradient, bool same_x)
    {
        KERNEL_TYPE::zero(gradient);
        jtf_func_->gra(&x[0],&gradient[0]);
    }

    virtual void eval_constraint_value(const vector_type &x, vector_type &c, bool same_x)
    {
        // ignore same_x
        size_t i = 0;
        value_type f = 0;
        for(; i < eqn_container_->size(); ++i) {
            f = 0;
            (*eqn_container_)[i]->val(&x[0], f);
            c[i] = f;
        }
        for(; i < eqn_container_->size() + ineqn_container_->size(); ++i) {
            f = 0;
            (*ineqn_container_)[i-eqn_container_->size()]->val(&x[0], f);
            c[i] = f;
        }
    }
    virtual void eval_constraint_gradient(const vector_type& x,
                                          typename KERNEL_TYPE::matrix_builder& c_gradient,
                                          bool same_x)
    {
        // ignore same x
        static std::vector<value_type> g;
        static std::vector<int_type> idx;
        KERNEL_TYPE::clear(c_gradient);
        size_t ci = 0, nnz;
        for(; ci < eqn_container_->size(); ++ci) {
            (*eqn_container_)[ci]->gra(&x[0], nnz, 0,0);
            g.resize(nnz);
            idx.resize(nnz);
            (*eqn_container_)[ci]->gra(&x[0], nnz, &g[0], &idx[0]);
            for(size_t i = 0; i < nnz; ++i) {
                c_gradient.push_back(std::make_pair(std::make_pair(ci, idx[i]),g[i]));
            }
        }
        for(; ci < eqn_container_->size() + ineqn_container_->size(); ++ci) {
            (*ineqn_container_)[ci-eqn_container_->size()]->gra(&x[0], nnz, 0,0);
            g.resize(nnz);
            idx.resize(nnz);
            (*ineqn_container_)[ci-eqn_container_->size()]->gra(&x[0], nnz, &g[0], &idx[0]);
            for(size_t i = 0; i < nnz; ++i) {
                c_gradient.push_back(std::make_pair(std::make_pair(ci, idx[i]),g[i]));
            }
        }
    }

private:
    void assemble(const vector_type &x)
    {
        {
            if(H_.size(1) == 0 || H_.size(2) == 0) {
                jtf_func_->hes(&x[0], nnz_h, format_h,0,0,0);
                if(nnz_h > 0) {
                    H_.resize(jtf_func_->dim(), jtf_func_->dim(), nnz_h);
                    jtf_func_->hes(&x[0], nnz_h, format_h, 0, &H_.ptr()[0], &H_.idx()[0]);
                }
            }

            if(!eqn_container_->empty()) {
                eqn_H_.resize(eqn_container_->size());
                size_t nnz , format;
                for(size_t eqi = 0; eqi < eqn_container_->size(); ++eqi) {
                    (*eqn_container_)[eqi]->hes(&x[0], nnz, format, 0,0,0);
                    if(nnz > 0) {
                        eqn_H_[eqi].resize(nr_inputs(),nr_inputs(), nnz);
                        (*eqn_container_)[eqi]->hes(&x[0],nnz,format, 0, &eqn_H_[eqi].ptr()[0],
                                                    &eqn_H_[eqi].idx()[0]);
                    }
                }
            }

            if(!ineqn_container_->empty()) {
                ineqn_H_.resize(ineqn_container_->size());
                size_t nnz , format;
                for(size_t eqi = 0; eqi < ineqn_container_->size(); ++eqi) {
                    (*ineqn_container_)[eqi]->hes(&x[0], nnz, format, 0,0,0);
                    if(nnz > 0) {
                        ineqn_H_[eqi].resize(nr_inputs(),nr_inputs(), nnz);
                        (*ineqn_container_)[eqi]->hes(&x[0],nnz,format, 0, &ineqn_H_[eqi].ptr()[0],
                                                      &ineqn_H_[eqi].idx()[0]);
                    }
                }
            }

            if(ij2double_ptr_hes_.empty()) {
                if(H_.size(1) != 0 || H_.size(2) != 0) {
                    for(size_t i = 0; i < H_.ptr().size() - 1; ++i) {
                        for(size_t pi = H_.ptr()[i]; pi < H_.ptr()[i+1]; ++pi) {
                            ij2double_ptr_hes_[std::make_pair(i,H_.idx()[pi])].push_back(
                                std::make_pair(-1,&H_.val()[pi]));
                        }
                    }
                }
                for(size_t eqi = 0; eqi < eqn_H_.size(); ++eqi) {
                    if(eqn_H_[eqi].size(1) == 0 || eqn_H_[eqi].size(2) == 0) continue;
                    for(size_t i = 0; i < eqn_H_[eqi].ptr().size() - 1; ++i) {
                        for(size_t pi = eqn_H_[eqi].ptr()[i]; pi < eqn_H_[eqi].ptr()[i+1]; ++pi) {
                            ij2double_ptr_hes_[std::make_pair(i,eqn_H_[eqi].idx()[pi])].push_back(
                                std::make_pair(eqi,&eqn_H_[eqi].val()[pi]));
                        }
                    }
                }
                for(size_t eqi = 0; eqi < ineqn_H_.size(); ++eqi) {
                    if(ineqn_H_[eqi].size(1) == 0 || ineqn_H_[eqi].size(2) == 0) continue;
                    for(size_t i = 0; i < ineqn_H_[eqi].ptr().size() - 1; ++i) {
                        for(size_t pi = ineqn_H_[eqi].ptr()[i]; pi < ineqn_H_[eqi].ptr()[i+1]; ++pi) {
                            ij2double_ptr_hes_[std::make_pair(i,ineqn_H_[eqi].idx()[pi])].push_back(
                                std::make_pair(eqi+eqn_container_->size(),&ineqn_H_[eqi].val()[pi]));
                        }
                    }
                }
            }
        }
    }
private:
    std::shared_ptr<JTF_FUNC> jtf_func_;
    std::shared_ptr<func_container> eqn_container_, ineqn_container_;
    hj::sparse::csc<value_type,int_type> H_;
    std::vector<hj::sparse::csc<value_type, int_type> > eqn_H_, ineqn_H_;
    std::map<std::pair<int_type,int_type>, std::vector<value_type*> >  ij2double_ptr_g_jac_;
    std::map<std::pair<int_type,int_type>, std::vector<std::pair<int_type,value_type*> > > ij2double_ptr_hes_;
    size_t nnz_h, format_h;
};

template <typename KERNEL_TYPE>
class jtf_func_to_obj_func_L_SR1 : public objective_function_L_SR1<KERNEL_TYPE>
{
public:
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;

    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef typename KERNEL_TYPE::int_type int_type;

    typedef jtf::function::functionN1_t<value_type, int_type> JTF_FUNC;
    typedef std::vector<std::shared_ptr<JTF_FUNC> > func_container;

    jtf_func_to_obj_func_L_SR1(std::shared_ptr<JTF_FUNC> jtf_obj_func,
                               std::shared_ptr<func_container> eqn_container,
                               std::shared_ptr<func_container> ineqn_constainer)
        :jtf_func_(jtf_obj_func), nnz_h(-1), format_h(-1) ,eqn_container_(eqn_container),
         ineqn_container_(ineqn_constainer)
    {
        if(!eqn_container_.get()) eqn_container_.reset(new func_container);
        if(!ineqn_container_.get()) ineqn_container_.reset(new func_container);
        static vector_type fake_x;
        KERNEL_TYPE::resize(nr_inputs(),fake_x);
        assemble(fake_x);
    }
    virtual ~jtf_func_to_obj_func_L_SR1() {}

    virtual size_t nr_inputs()const
    {
        return jtf_func_->dim();
    }

    virtual size_t nr_equals()const
    {
        if(!eqn_container_.get()) return 0;
        return eqn_container_->size();
    }

    virtual size_t nr_inequals()const
    {
        if(!ineqn_container_.get()) return 0;
        return ineqn_container_->size();
    }

    virtual value_type eval(const vector_type& x,bool same_x)
    {
        // ignore same_x
        value_type f = 0;
        jtf_func_->val(&x[0], f);
        return f;
    }

    virtual void eval_func_gradient(const vector_type &x, vector_type &gradient, bool same_x)
    {
        KERNEL_TYPE::zero(gradient);
        jtf_func_->gra(&x[0],&gradient[0]);
    }

    virtual void eval_lagrangian_hessian(const vector_type& x,
                                         const vector_type& lambda,
                                         value_type coef,
                                         typename KERNEL_TYPE::matrix_builder& hessian,
                                         bool same_x)
    {
        // ignore same_x
        size_t format = 1;
        H_.val() *= 0;
        jtf_func_->hes(&x[0], nnz_h, format_h, &H_.val()[0], &H_.ptr()[0], &H_.idx()[0]);
        H_.val() *= coef;

        for(size_t eqi = 0; eqi < eqn_container_->size(); ++eqi) {
            if(eqn_H_[eqi].size(1) != 0 && eqn_H_[eqi].size(2) != -1) {
                eqn_H_[eqi].val() *= 0;
                size_t nnz = hj::sparse::nnz(eqn_H_[eqi]);
                if(nnz > 0)
                    (*eqn_container_)[eqi]->hes(&x[0], nnz , format, 0,
                                                &eqn_H_[eqi].ptr()[0], &eqn_H_[eqi].idx()[0]);
            }
        }

        for(size_t eqi = 0; eqi < ineqn_container_->size(); ++eqi) {
            if(ineqn_H_[eqi].size(1) != 0 && ineqn_H_[eqi].size(2) != -1) {
                ineqn_H_[eqi].val() *= 0;
                size_t nnz = hj::sparse::nnz(ineqn_H_[eqi]);
                if(nnz > 0)
                    (*ineqn_container_)[eqi]->hes(&x[0], nnz , format, &ineqn_H_[eqi].val()[0],
                                                  &ineqn_H_[eqi].ptr()[0], &ineqn_H_[eqi].idx()[0]);
            }
        }

        KERNEL_TYPE::clear(hessian);
        if(!ij2double_ptr_hes_.empty()) {
            //for(const auto & term : ij2double_ptr_hes_) {
            for(typename std::map<std::pair<int_type,int_type>, std::vector<std::pair<int_type,value_type*> > >::const_iterator
                    cit = ij2double_ptr_hes_.begin(); cit != ij2double_ptr_hes_.end(); ++cit) {
                const std::pair<std::pair<int_type,int_type>, std::vector<std::pair<int_type,value_type*> > > & term = *cit;
                value_type h_temp = 0;
                // for(const auto & ptr : term.second) {
                for(typename std::vector<std::pair<int_type,value_type*> >::const_iterator pit = term.second.begin();
                        pit != term.second.end(); ++pit) {
                    const std::pair<int_type,value_type*>  & ptr = *pit;
                    double v = *ptr.second;
                    if(ptr.first != -1)// if is constraint
                        v *= lambda[ptr.first];
                    h_temp += v;
                }

                hessian.push_back(std::make_pair(std::make_pair(term.first.first, term.first.second), h_temp));
            }
        }
    }
    virtual void eval_constraint_value(const vector_type &x, vector_type &c, bool same_x)
    {
        // ignore same_x
        size_t i = 0;
        value_type f = 0;
        for(; i < eqn_container_->size(); ++i) {
            f = 0;
            (*eqn_container_)[i]->val(&x[0], f);
            c[i] = f;
        }
        for(; i < eqn_container_->size() + ineqn_container_->size(); ++i) {
            f = 0;
            (*ineqn_container_)[i-eqn_container_->size()]->val(&x[0], f);
            c[i] = f;
        }
    }
    virtual void eval_constraint_gradient(const vector_type& x,
                                          typename KERNEL_TYPE::matrix_builder& c_gradient,
                                          bool same_x)
    {
        // ignore same x
        static std::vector<value_type> g;
        static std::vector<int_type> idx;
        KERNEL_TYPE::clear(c_gradient);
        size_t ci = 0, nnz;
        for(; ci < eqn_container_->size(); ++ci) {
            (*eqn_container_)[ci]->gra(&x[0], nnz, 0,0);
            g.resize(nnz);
            idx.resize(nnz);
            (*eqn_container_)[ci]->gra(&x[0], nnz, &g[0], &idx[0]);
            for(size_t i = 0; i < nnz; ++i) {
                c_gradient.push_back(std::make_pair(std::make_pair(ci, idx[i]),g[i]));
            }
        }
        for(; ci < eqn_container_->size() + ineqn_container_->size(); ++ci) {
            (*ineqn_container_)[ci-eqn_container_->size()]->gra(&x[0], nnz, 0,0);
            g.resize(nnz);
            idx.resize(nnz);
            (*ineqn_container_)[ci-eqn_container_->size()]->gra(&x[0], nnz, &g[0], &idx[0]);
            for(size_t i = 0; i < nnz; ++i) {
                c_gradient.push_back(std::make_pair(std::make_pair(ci, idx[i]),g[i]));
            }
        }
    }

private:
    void assemble(const vector_type &x)
    {
        {
            if(H_.size(1) == 0 || H_.size(2) == 0) {
                jtf_func_->hes(&x[0], nnz_h, format_h,0,0,0);
                if(nnz_h > 0) {
                    H_.resize(jtf_func_->dim(), jtf_func_->dim(), nnz_h);
                    jtf_func_->hes(&x[0], nnz_h, format_h, 0, &H_.ptr()[0], &H_.idx()[0]);
                }
            }

            if(!eqn_container_->empty()) {
                eqn_H_.resize(eqn_container_->size());
                size_t nnz , format;
                for(size_t eqi = 0; eqi < eqn_container_->size(); ++eqi) {
                    (*eqn_container_)[eqi]->hes(&x[0], nnz, format, 0,0,0);
                    if(nnz > 0) {
                        eqn_H_[eqi].resize(nr_inputs(),nr_inputs(), nnz);
                        (*eqn_container_)[eqi]->hes(&x[0],nnz,format, 0, &eqn_H_[eqi].ptr()[0],
                                                    &eqn_H_[eqi].idx()[0]);
                    }
                }
            }

            if(!ineqn_container_->empty()) {
                ineqn_H_.resize(ineqn_container_->size());
                size_t nnz , format;
                for(size_t eqi = 0; eqi < ineqn_container_->size(); ++eqi) {
                    (*ineqn_container_)[eqi]->hes(&x[0], nnz, format, 0,0,0);
                    if(nnz > 0) {
                        ineqn_H_[eqi].resize(nr_inputs(),nr_inputs(), nnz);
                        (*ineqn_container_)[eqi]->hes(&x[0],nnz,format, 0, &ineqn_H_[eqi].ptr()[0],
                                                      &ineqn_H_[eqi].idx()[0]);
                    }
                }
            }

            if(ij2double_ptr_hes_.empty()) {
                if(H_.size(1) != 0 || H_.size(2) != 0) {
                    for(size_t i = 0; i < H_.ptr().size() - 1; ++i) {
                        for(size_t pi = H_.ptr()[i]; pi < H_.ptr()[i+1]; ++pi) {
                            ij2double_ptr_hes_[std::make_pair(i,H_.idx()[pi])].push_back(
                                std::make_pair(-1,&H_.val()[pi]));
                        }
                    }
                }
                for(size_t eqi = 0; eqi < eqn_H_.size(); ++eqi) {
                    if(eqn_H_[eqi].size(1) == 0 || eqn_H_[eqi].size(2) == 0) continue;
                    for(size_t i = 0; i < eqn_H_[eqi].ptr().size() - 1; ++i) {
                        for(size_t pi = eqn_H_[eqi].ptr()[i]; pi < eqn_H_[eqi].ptr()[i+1]; ++pi) {
                            ij2double_ptr_hes_[std::make_pair(i,eqn_H_[eqi].idx()[pi])].push_back(
                                std::make_pair(eqi,&eqn_H_[eqi].val()[pi]));
                        }
                    }
                }
                for(size_t eqi = 0; eqi < ineqn_H_.size(); ++eqi) {
                    if(ineqn_H_[eqi].size(1) == 0 || ineqn_H_[eqi].size(2) == 0) continue;
                    for(size_t i = 0; i < ineqn_H_[eqi].ptr().size() - 1; ++i) {
                        for(size_t pi = ineqn_H_[eqi].ptr()[i]; pi < ineqn_H_[eqi].ptr()[i+1]; ++pi) {
                            ij2double_ptr_hes_[std::make_pair(i,ineqn_H_[eqi].idx()[pi])].push_back(
                                std::make_pair(eqi+eqn_container_->size(),&ineqn_H_[eqi].val()[pi]));
                        }
                    }
                }
            }
        }
    }
private:
    std::shared_ptr<JTF_FUNC> jtf_func_;
    std::shared_ptr< func_container> eqn_container_, ineqn_container_;
    hj::sparse::csc<value_type,int_type> H_;
    std::vector<hj::sparse::csc<value_type, int_type> > eqn_H_, ineqn_H_;
    std::map<std::pair<int_type,int_type>, std::vector<value_type*> >  ij2double_ptr_g_jac_;
    std::map<std::pair<int_type,int_type>, std::vector<std::pair<int_type,value_type*> > > ij2double_ptr_hes_;
    size_t nnz_h, format_h;
};
}
}

#endif // JTF_FUNCTION_TO_OBJ_FUNC_H
