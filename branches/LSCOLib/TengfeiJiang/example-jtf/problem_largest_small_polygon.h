#ifndef PROBLEM_LARGEST_SMALL_POLYGON_H
#define PROBLEM_LARGEST_SMALL_POLYGON_H

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>

#include <jtflib/function/function.h>
#include <jtflib/function/func_aux.h>
#include "jtf_func_to_obj_func.h"

#include <iostream>
#include <fstream>
#include "kernel_zjucad_matrix.h"
#include "derivative_checker.h"
//#include "third_solver_bleic.h"
#include "third_solver_ipopt.h"

namespace zjucad
{
namespace LSCO
{
template <typename MATRIX_TYPE>
struct kernel_traits;

typedef jtf::function::functionN1_t<double, int32_t> jtf_func;
typedef std::shared_ptr<jtf_func> jtf_func_cons_ptr;
typedef std::vector<jtf_func_cons_ptr> func_container;

class largest_small_polygon_func : public jtf_func
{
    // y = 0.5*r_i+1*r_i*sin(theta_i-theta_i+1)
public:
    typedef jtf_func::val_type val_type;
    typedef jtf_func::int_type int_type;

    largest_small_polygon_func(const size_t variable_num,
                               const size_t idx)
        :vn_(variable_num), idx_(idx)
    {
        assert(idx_ *2 < vn_);
    }
    virtual ~largest_small_polygon_func() {}
    virtual size_t dim() const
    {
        return vn_;
    }
    virtual int val(const val_type *x, val_type &v) const
    {
        using namespace zjucad::matrix;
        itr_matrix<const val_type*> x0(dim()/2,2,x); // x = [r, theta]
        v += 0.5*(x0(idx_,0)*x0(idx_+1,0)*std::sin(x0(idx_,1) - x0(idx_+1,1)));
        return 0;
    }
    virtual int gra(const val_type *x, size_t &nnz, val_type * g, int_type *idx)
    {
        using namespace zjucad::matrix;
        if(g == 0 && idx == 0) {
            nnz = 4;
            return 0;
        }
        itr_matrix<const val_type*> x0(dim()/2,2,x); // x = [r, theta]
        g[0] = 0.5*x0(idx_+1,0)*std::sin(x0(idx_,1)-x0(idx_+1,1));
        idx[0] = idx_;

        g[1] = 0.5*x0(idx_,0)*std::sin(x0(idx_,1)-x0(idx_+1,1));
        idx[1] = idx_+1;

        g[2] = 0.5*x0(idx_+1,0)*x0(idx_,0)*std::cos(x0(idx_,1)-x0(idx_+1,1));
        idx[2] = idx_ + x0.size(1);

        g[3] = -0.5*x0(idx_+1,0)*x0(idx_,0)*std::cos(x0(idx_,1)-x0(idx_+1,1));
        idx[3] = idx_+1+x0.size(1);

        return 0;
    }
    virtual int gra(const val_type *x, val_type *g)
    {
        static val_type g_[4];
        static int_type idx_[4];
        size_t nnz;
        gra(x, nnz, &g_[0], &idx_[0]);
        for(size_t i = 0; i < 4; ++i) {
            g[idx_[i]] += g_[i];
        }
        return 0;
    }
    virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                    val_type *h, int_type *ptr, int_type *idx, val_type alpha = 1)
    {
        if(h == 0 && idx == 0 && ptr == 0) {
            nnz = 14;
            format = 1;
            return 0;
        }
        const int_type theta_i_idx = vn_/2 + idx_;
        const int_type theta_i_1_idx = vn_/2 + idx_+1;
        const int_type r_i_idx = idx_;
        const int_type r_i_idx_1 = idx_+1;
        const val_type & theta_i = x[theta_i_idx];
        const val_type & theta_i_1 = x[theta_i_1_idx];
        const val_type & r_i = x[r_i_idx];
        const val_type & r_i_1 = x[r_i_idx_1];

        if(h == 0 && idx != 0 && ptr != 0) {
            ptr[r_i_idx_1] = ptr[r_i_idx]+ 3;
            idx[ptr[r_i_idx]+0] = r_i_idx_1;
            idx[ptr[r_i_idx]+1] = theta_i_idx;
            idx[ptr[r_i_idx]+2] = theta_i_1_idx;

            ptr[r_i_idx_1+1] = ptr[r_i_idx_1]+ 3;
            idx[ptr[r_i_idx_1]+0] = r_i_idx;
            idx[ptr[r_i_idx_1]+1] = theta_i_idx;
            idx[ptr[r_i_idx_1]+2] = theta_i_1_idx;

            ptr[theta_i_idx] = ptr[r_i_idx_1+1];

            ptr[theta_i_1_idx] = ptr[theta_i_idx]+ 4;
            idx[ptr[theta_i_idx]+0] = r_i_idx;
            idx[ptr[theta_i_idx]+1] = r_i_idx_1;
            idx[ptr[theta_i_idx]+2] = theta_i_idx;
            idx[ptr[theta_i_idx]+3] = theta_i_1_idx;

            ptr[theta_i_1_idx+1] = ptr[theta_i_1_idx]+ 4;
            idx[ptr[theta_i_1_idx]+0] = r_i_idx;
            idx[ptr[theta_i_1_idx]+1] = r_i_idx_1;
            idx[ptr[theta_i_1_idx]+2] = theta_i_idx;
            idx[ptr[theta_i_1_idx]+3] = theta_i_1_idx;
            return 0;
        }

        if(h != 0 && idx != 0 && ptr != 0) {
            const val_type sin_theta_i_theta_i_1 = std::sin(theta_i-theta_i_1);
            const val_type cos_theta_i_theta_i_1 = std::cos(theta_i-theta_i_1);

            jtf::function::add_to_csc(h,ptr,idx,r_i_idx, r_i_idx_1, 0.5*sin_theta_i_theta_i_1);
            jtf::function::add_to_csc(h,ptr,idx,r_i_idx, theta_i_idx, 0.5*r_i_1* cos_theta_i_theta_i_1);
            jtf::function::add_to_csc(h,ptr,idx,r_i_idx, theta_i_1_idx, -0.5*r_i_1* cos_theta_i_theta_i_1);

            jtf::function::add_to_csc(h,ptr,idx,r_i_idx_1, r_i_idx, 0.5*sin_theta_i_theta_i_1);
            jtf::function::add_to_csc(h,ptr,idx,r_i_idx_1, theta_i_idx, 0.5*r_i* cos_theta_i_theta_i_1);
            jtf::function::add_to_csc(h,ptr,idx,r_i_idx_1, theta_i_1_idx, -0.5*r_i* cos_theta_i_theta_i_1);

            jtf::function::add_to_csc(h,ptr,idx,theta_i_idx, r_i_idx, 0.5*r_i_1*cos_theta_i_theta_i_1);
            jtf::function::add_to_csc(h,ptr,idx,theta_i_idx, r_i_idx_1, 0.5*r_i*cos_theta_i_theta_i_1);
            jtf::function::add_to_csc(h,ptr,idx,theta_i_idx, theta_i_idx, -0.5*r_i*r_i_1* sin_theta_i_theta_i_1);
            jtf::function::add_to_csc(h,ptr,idx,theta_i_idx, theta_i_1_idx, -0.5*r_i_1* cos_theta_i_theta_i_1);

            jtf::function::add_to_csc(h,ptr,idx,theta_i_1_idx, r_i_idx, -0.5*r_i_1*cos_theta_i_theta_i_1);
            jtf::function::add_to_csc(h,ptr,idx,theta_i_1_idx, r_i_idx_1, -0.5*r_i*cos_theta_i_theta_i_1);
            jtf::function::add_to_csc(h,ptr,idx,theta_i_1_idx, theta_i_idx, 0.5*r_i*r_i_1* sin_theta_i_theta_i_1);
            jtf::function::add_to_csc(h,ptr,idx,theta_i_1_idx, theta_i_1_idx, -0.5*r_i*r_i_1* sin_theta_i_theta_i_1);
        }

        return 0;
    }
    virtual int hes_block(const val_type *x, val_type *h, val_type alpha)
    {
        return __LINE__;
    }
private:
    const size_t vn_;
    const size_t idx_;
};

class first_constraint_func : public jtf_func
{
    // 1 - r_i^2-r_j^2+2*r_i*r_j*std::cos(theta_i-theta_j)
public:
    first_constraint_func(const size_t variable_number, const size_t i, const size_t j)
        :vn_(variable_number), i_(i), j_(j)
    {
        assert(i < j);
    }
    virtual size_t dim() const
    {
        return vn_;
    }
    virtual int val(const val_type *x, val_type &v) const
    {
        using namespace zjucad::matrix;
        itr_matrix<const val_type*> x0(vn_/2,2, x);
        const val_type & ri = x0(i_,0);
        const val_type & rj = x0(j_,0);
        const val_type & thetai = x0(i_,1);
        const val_type & thetaj = x0(j_,1);
        v += 1-ri*ri-rj*rj+2*ri*rj*std::cos(thetai-thetaj);
        return 0;
    }
    virtual int gra(const val_type *x, val_type *g)
    {
        static val_type g_[4];
        static int_type idx_[4];
        size_t nnz = 0;
        gra(x,nnz, g_, idx_);
        for(size_t i = 0; i < 4; ++i) {
            g[idx_[i]] += g_[i];
        }
        return 0;
    }
    virtual int gra(const val_type *x, size_t &nnz, val_type *g, int_type *idx)
    {
        if(g == 0 && idx == 0) {
            nnz = 4;
            return 0;
        }
        using namespace zjucad::matrix;
        itr_matrix<const val_type*> x0(vn_/2,2, x);
        const val_type & ri = x0(i_,0);
        const val_type & rj = x0(j_,0);
        const val_type & thetai = x0(i_,1);
        const val_type & thetaj = x0(j_,1);
        const val_type cos_thetai_thetaj = std::cos(thetai-thetaj);
        const val_type sin_thetai_thetaj = std::sin(thetai - thetaj);

        g[0] = -1*(2*ri-2*rj*cos_thetai_thetaj);
        idx[0] = i_;

        g[1] = -1*(2*rj-2*ri*cos_thetai_thetaj);
        idx[1] = j_;

        g[2] = -2*ri*rj*sin_thetai_thetaj;
        idx[2] = i_ + vn_/2;

        g[3] = 2*ri*rj*sin_thetai_thetaj;
        idx[3] = j_ + vn_/2;

        return 0;
    }
    virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                    val_type *h, int_type *ptr, int_type *idx, val_type alpha=1)
    {
        nnz = 0;
        return 0;
    }
    virtual int hes_block(const val_type *x, val_type *h, val_type alpha)
    {
        return __LINE__;
    }
private:
    const size_t i_, j_;
    size_t vn_;
};

class theta_constraint_0_func : public jtf_func
{
    // theta_i+1 -theta_i >= 0
public:
    theta_constraint_0_func(const size_t variable_number,
                            const size_t i)
        :vn_(variable_number), i_(i) {}
    virtual ~theta_constraint_0_func() {}
    virtual size_t dim()const
    {
        return vn_;
    }
    virtual int val(const val_type *x, val_type &v) const
    {
        v += x[i_+vn_/2+1] - x[i_+vn_/2];
        return 0;
    }
    virtual int gra(const val_type * x, val_type *g)
    {
        static val_type g_[2];
        static int_type idx_[2];
        size_t nnz;
        gra(x, nnz, g_,idx_);
        for(size_t i = 0; i < 2; ++i)
            g[idx_[i]] == g_[i];
        return 0;
    }
    virtual int gra(const val_type *x, size_t &nnz, val_type *g, int_type *idx)
    {
        if(g == 0 && idx == 0) {
            nnz = 2;
            return 0;
        }
        g[0] = 1;
        idx[0] = i_+vn_/2+1;

        g[1] = -1;
        idx[1] = i_+vn_/2;
        return 0;
    }
    virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                    val_type *h, int_type *ptr, int_type *idx, val_type alpha)
    {
        nnz = 0;
        return __LINE__;
    }
    virtual int hes_block(const val_type *x, val_type *h, val_type alpha)
    {
        return __LINE__;
    }
private:
    const size_t vn_;
    const size_t i_;
};

class theta_bound_constraint_0_func : public jtf_func
{
    // theta_i >= 0
public:
    theta_bound_constraint_0_func(const size_t variable_number,
                                  const size_t i)
        : vn_(variable_number), i_(i) {}
    virtual ~theta_bound_constraint_0_func() {}
    virtual size_t dim() const
    {
        return vn_;
    }
    virtual int val(const val_type *x, val_type &v) const
    {
        v += x[i_+vn_/2];
        return 0;
    }
    virtual int gra(const val_type *x, val_type *g)
    {
        g[i_+vn_/2] += 1.0;
        return 0;
    }
    virtual int gra(const val_type *x, size_t & nnz, val_type *g, int_type *idx)
    {
        if(g == 0 && idx == 0) {
            nnz = 1;
            return 0;
        }
        g[0] = 1;
        idx[0] = i_+vn_/2;
        return 0;
    }
    virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                    val_type *h, int_type *ptr, int_type *idx, val_type alpha)
    {
        nnz = 0;
        return __LINE__;
    }
    virtual int hes_block(const val_type *x, val_type *h, val_type alpha)
    {
        return __LINE__;
    }
private:
    const size_t vn_;
    const size_t i_;
};

class theta_bound_constraint_1_func : public jtf_func
{
    // pi-theta_i
public:
    theta_bound_constraint_1_func(const size_t variable_number,
                                  const size_t i)
        : vn_(variable_number), i_(i) {}
    virtual ~theta_bound_constraint_1_func() {}
    virtual size_t dim() const
    {
        return vn_;
    }
    virtual int val(const val_type *x, val_type &v) const
    {
        v += 3.1415925535897-x[i_+vn_/2];
        return 0;
    }
    virtual int gra(const val_type *x, val_type *g)
    {
        g[i_+vn_/2] += -1.0;
        return 0;
    }
    virtual int gra(const val_type *x, size_t & nnz, val_type *g, int_type *idx)
    {
        if(g == 0 && idx == 0) {
            nnz = 1;
            return 0;
        }
        g[0] = -1;
        idx[0] = i_+vn_/2;
        return 0;
    }
    virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                    val_type *h, int_type *ptr, int_type *idx, val_type alpha)
    {
        nnz = 0;
        return __LINE__;
    }
    virtual int hes_block(const val_type *x, val_type *h, val_type alpha)
    {
        return __LINE__;
    }
private:
    const size_t vn_;
    const size_t i_;
};

class r_bound_constraint_func : public jtf_func
{
    // ri>= 0
public:
    r_bound_constraint_func(const size_t variable_number,
                            const size_t idx):vn_(variable_number),idx_(idx) {}
    virtual ~r_bound_constraint_func() {}
    virtual size_t dim() const
    {
        return vn_;
    }
    virtual int val(const val_type *x, val_type &v) const
    {
        v += x[idx_];
        return 0;
    }
    virtual int gra(const val_type *x, val_type *g)
    {
        g[idx_] += 1.0;
        return 0;
    }
    virtual int gra(const val_type *x, size_t &nnz, val_type*g, int_type*idx)
    {
        if(g == 0 && idx == 0) {
            nnz = 1;
            return 0;
        }
        g[0] = 1.0;
        idx[0] = idx_;
        return 0;
    }
    virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                    val_type *h, int_type *ptr, int_type *idx, val_type alpha)
    {
        nnz = 0;
        return __LINE__;
    }
    virtual int hes_block(const val_type *x, val_type *h, val_type alpha)
    {
        return __LINE__;
    }
private:
    const size_t vn_,idx_;
};

class r_bound_constraint_1_func : public jtf_func
{
    // 1-ri >= 0
public:
    r_bound_constraint_1_func(const size_t variable_number,
                              const size_t idx):vn_(variable_number),idx_(idx) {}
    virtual ~r_bound_constraint_1_func() {}
    virtual size_t dim() const
    {
        return vn_;
    }
    virtual int val(const val_type *x, val_type &v) const
    {
        v += 1-x[idx_];
        return 0;
    }
    virtual int gra(const val_type *x, val_type *g)
    {
        g[idx_] += -1.0;
        return 0;
    }
    virtual int gra(const val_type *x, size_t &nnz, val_type*g, int_type*idx)
    {
        if(g == 0 && idx == 0) {
            nnz = 1;
            return 0;
        }
        g[0] = -1.0;
        idx[0] = idx_;
        return 0;
    }
    virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                    val_type *h, int_type *ptr, int_type *idx, val_type alpha)
    {
        nnz = 0;
        return __LINE__;
    }
    virtual int hes_block(const val_type *x, val_type *h, val_type alpha)
    {
        return __LINE__;
    }
private:
    const size_t vn_,idx_;
};


template <typename OS, typename FLOAT, typename INT>
void line2vtk(
    OS &os,
    const FLOAT *node, size_t node_num,
    const INT *line, size_t line_num)
{
    os << "# vtk DataFile Version 2.0\nTRI\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";

    os<< "POINTS " << node_num << " float\n";
    for(size_t i = 0; i < node_num; ++i)
        os << node[i*3+0] << " " << node[i*3+1] << " " << node[i*3+2] << "\n";

    os << "CELLS " << line_num << " " << line_num*3 << "\n";
    for(size_t i = 0; i < line_num; ++i)
        os << 2 << " " << line[i*2+0] << " " << line[i*2+1] << "\n";

    os << "CELL_TYPES " << line_num << "\n";
    for(size_t i = 0; i < line_num; ++i)
        os << 3 << "\n";
}


void to_vtk(const zjucad::matrix::matrix<double> & r_theta,
            const char * output_file)
{
    using namespace zjucad::matrix;
    zjucad::matrix::matrix<double> points = zeros<double>(3, r_theta.size(1));
    for(size_t pi = 0; pi < points.size(2); ++pi) {
        points(0,pi) = r_theta(pi,0) * std::cos(r_theta(pi,1));
        points(1,pi) = r_theta(pi,0) * std::sin(r_theta(pi,1));
    }

    std::vector<size_t> lines(r_theta.size(1)*2);
    for(size_t i = 0; i < r_theta.size(1); ++i) {
        lines[i*2+0] = i;
        lines[i*2+1] = (i+1)%r_theta.size(1);
    }

    std::ofstream ofs(output_file);
    line2vtk(ofs, &points[0], points.size(2), &lines[0], lines.size()/2);
}

void init_r_theta(zjucad::matrix::matrix<double> & r_theta,
                  const size_t n_side)
{
    const double angle = 2*3.1415926535897/n_side;
    double coord[2] = {0,0};
    for(size_t i = 0; i < n_side; ++i) {
        coord[0] = 0.5 * sin((i+1)*angle);
        coord[1] = 0.5-0.5*cos((i+1)*angle);
        r_theta(i,0) = std::sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
        r_theta(i,1) = std::atan2(coord[1],coord[0]);
    }
}

void largest_small_polygon(const size_t n_side)
{
    typedef std::vector<jtf_func_cons_ptr> jtf_func_container;

    zjucad::matrix::matrix<double> r_theta = zjucad::matrix::zeros<double>(n_side,2);
    init_r_theta(r_theta, n_side);
    to_vtk(r_theta, "init_largest_area.vtk");
    std::cerr << r_theta(n_side-1,zjucad::matrix::colon()) << std::endl;

    jtf_func_container jfc;
    {
        for(size_t i = 0; i < n_side-1; ++i) {
            jfc.push_back(jtf_func_cons_ptr(new largest_small_polygon_func(2*n_side, i)));
        }
    }

    jtf_func_cons_ptr obj(new jtf::function::sum_function<double,int32_t>(jfc));

    std::shared_ptr<func_container > eqn_func(new func_container),
        ineqn_func(new func_container),test(new func_container);

    {
        for(size_t i = 0; i < n_side; ++i) {
            for(size_t j = i + 1; j < n_side; ++j) {
                ineqn_func->push_back(jtf_func_cons_ptr(new first_constraint_func(2*n_side,i,j)));
            }
        }

        for(size_t i = 0; i < n_side-1; ++i) {
            ineqn_func->push_back(jtf_func_cons_ptr(new theta_constraint_0_func(2*n_side,i)));
        }

        for(size_t i = 0; i < n_side; ++i) {
            ineqn_func->push_back(jtf_func_cons_ptr(new theta_bound_constraint_0_func(2*n_side,i)));
            ineqn_func->push_back(jtf_func_cons_ptr(new theta_bound_constraint_1_func(2*n_side,i)));
            ineqn_func->push_back(jtf_func_cons_ptr(new r_bound_constraint_func(2*n_side,i)));
        }

        for(size_t i = 0; i < n_side-1; ++i) {
            ineqn_func->push_back(jtf_func_cons_ptr(new r_bound_constraint_1_func(2*n_side,i)));
        }

    }

    {
        eqn_func->push_back(jtf_func_cons_ptr(new theta_bound_constraint_1_func(2*n_side,n_side-1)));
        eqn_func->push_back(jtf_func_cons_ptr(new r_bound_constraint_func(2*n_side, n_side-1)));
    }

    typedef kernel_traits<typename hj::sparse::csc<double,int32_t> > kernel_type;

    boost::shared_ptr<objective_function<kernel_type> > jftof(
        new jtf_func_to_obj_func<kernel_type>(obj,eqn_func, ineqn_func));

    //      {
    //        double total_area = 0;
    //        obj->val(&r_theta[0], total_area);
    //        std::cerr << "init area = " << -1*total_area << std::endl;

    //        std::cerr << "# check equality constraints " << std::endl;
    //        double v = 0;
    //        for(size_t ei = 0; ei < eqn_func->size(); ++ei){
    //            v = 0;
    //            (*eqn_func)[ei]->val(&r_theta[0], v);
    //            std::cerr << "+ "  << ei << " " << v << std::endl;
    //          }

    //        std::cerr << "# check inequality constraints " << std::endl;
    //        v = 0;
    //        for(size_t ei = 0; ei < ineqn_func->size(); ++ei){
    //            v = 0;
    //            (*ineqn_func)[ei]->val(&r_theta[0], v);
    //            std::cerr << "+ "  << ei << " " << v << std::endl;
    //          }
    //      }

    third_solver_ipopt tsb(jftof);
    boost::property_tree::ptree pt;
    pt.put("iter.value",2000);
    tsb.solve(r_theta, pt);

    to_vtk(r_theta, "largest_area.vtk");
    double total_area = 0;
    obj->val(&r_theta[0], total_area);
    std::cerr << "max area = " << -1*total_area << std::endl;
    //std::cerr << r_theta(zjucad::matrix::colon(3,5),zjucad::matrix::colon()) << std::endl;

    //      {
    //        double total_area = 0;
    //        obj->val(&r_theta[0], total_area);
    //        std::cerr << "max area = " << -1*total_area << std::endl;

    //        std::cerr << "# check equality constraints " << std::endl;
    //        double v = 0;
    //        for(size_t ei = 0; ei < eqn_func->size(); ++ei){
    //            v = 0;
    //            (*eqn_func)[ei]->val(&r_theta[0], v);
    //            std::cerr << "+ "  << ei << " " << v << std::endl;
    //          }

    //        std::cerr << "# check inequality constraints " << std::endl;
    //        v = 0;
    //        for(size_t ei = 0; ei < ineqn_func->size(); ++ei){
    //            v = 0;
    //            (*ineqn_func)[ei]->val(&r_theta[0], v);
    //            std::cerr << "+ "  << ei << " " << v << std::endl;
    //          }
    //      }
}
}
}
#endif // PROBLEM_DEBUG_H
