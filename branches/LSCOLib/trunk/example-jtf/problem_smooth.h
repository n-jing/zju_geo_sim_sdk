#ifndef PROBLEM_SMOOTH_H
#define PROBLEM_SMOOTH_H

#include <fstream>
#include <jtflib/function/function.h>
#include <jtflib/function/func_aux.h>

#include "third_solver_ipopt.h"
#include "jtf_func_to_obj_func.h"
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>

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
typedef std::shared_ptr<hj_func> hj_func_cons_ptr;

class node_fix2 : public jtf_func
{
public :
    node_fix2(const size_t node_number,
              const size_t idx,
              const zjucad::matrix::matrix<val_type> & fix_point,
              const double w)
        :node_number_(node_number), idx_(idx), fix_point_(fix_point), w_(w) {}
    virtual ~node_fix2() {}
    virtual size_t dim() const {
        return 3 * node_number_;
    }
    virtual int val(const val_type *x, val_type &v) const {
        using namespace zjucad::matrix;
        zjucad::matrix::itr_matrix<const val_type *> x0(3, node_number_,x);
        zjucad::matrix::matrix<val_type> diff = x0(colon(),idx_) - fix_point_;
        v += w_*0.5*dot(diff,diff);
        return 0;
    }
    virtual int gra(const val_type *x, size_t &nnz, val_type *g, int_type *idx) {
        if(g == 0 && idx == 0) {
            nnz = 3;
            return 0;
        }
        using namespace zjucad::matrix;
        zjucad::matrix::itr_matrix<const val_type *> x0(3, node_number_,x);
        zjucad::matrix::itr_matrix<val_type*> g_m(3,1,g);
        g_m = x0(colon(),idx_) - fix_point_;
        g_m *= w_;
        for(size_t i = 0; i < 3; ++i)
            idx[i] = 3*idx_+i;

        return 0;
    }
    virtual int gra(const val_type *x, val_type *g) {
        using namespace zjucad::matrix;
        zjucad::matrix::itr_matrix<const val_type *> x0(3, node_number_,x);
        for(size_t i = 0 ; i < 3; ++i) {
            g[3*idx_+i] += w_*(x0(i,idx_) - fix_point_[i]);
        }
        return 0;
    }

    virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                    val_type *h, int_type *ptr, int_type *idx, val_type alpha) {
        if(h == 0  && idx == 0 && ptr == 0) {
            nnz = 3;
            format = 1;
            return 0;
        }
        if(h == 0 && idx != 0 && ptr != 0) {
            for(size_t di = 0; di < 3; ++di) {
                ptr[3*idx_+di+1] = ptr[3*idx_+di] + 1;
                idx[ptr[3*idx_+di]]  = 3*idx_+di;
            }
            return 0;
        }
        if(h != 0 && idx != 0 && ptr != 0) {
            for(size_t di = 0; di < 3; ++di)
                jtf::function::add_to_csc(h,ptr,idx, 3*idx_+di, 3*idx_+di,w_);
        }
        return 0;
    }
    virtual int hes_block(const val_type *x, val_type *h, val_type alpha) {
        return __LINE__;
    }
private:
    const size_t node_number_;
    const size_t idx_;
    const zjucad::matrix::matrix<val_type> fix_point_;
    const double w_;
};

class degree_smooth_func : public hj_func
{
public:
    degree_smooth_func(const size_t point_num,
                       const size_t point_idx,
                       const std::vector<size_t> & one_ring)
        :point_num_(point_num), point_idx_(point_idx), one_ring_(one_ring), dim_(3) {}
    virtual ~degree_smooth_func() {}

public:
    virtual size_t dim_of_x(void) const
    {
        return point_num_ * dim_;
    }
    virtual size_t dim_of_f(void) const
    {
        return dim_;
    }
    virtual int val(const double *x, double *f,
                    hj::function::func_ctx *ctx = 0)const
    {
        assert(x);
        using namespace zjucad::matrix;
        zjucad::matrix::itr_matrix<double* > f_m(dim_,1,f);
        zjucad::matrix::matrix<double> one_ring_node;
        zjucad::matrix::itr_matrix<const double *> x_m(dim_, point_num_, x);
        itr_matrix<const size_t*> one_ring_mat(one_ring_.size(),1, &one_ring_[0]);
        one_ring_node = x_m(zjucad::matrix::colon(), one_ring_mat);

        f_m  *= 0;
        for(size_t i = 0; i < one_ring_.size(); ++i) {
            f_m -= one_ring_node(colon(),i);
        }
        zjucad::matrix::matrix<double> point_(dim_,1);

        point_ = x_m(zjucad::matrix::colon(), point_idx_);

        f_m /= one_ring_.size();
        f_m += point_;
        return 0;
    }
    virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                    int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const
    {

        for(size_t di = 0; di < dim_; ++di) {
            ptr[di+1] = ptr[di] + one_ring_.size() + 1;
            for(size_t i = 0; i < one_ring_.size(); ++i) {
                idx[ptr[di] + i] = dim_*one_ring_[i]+di;
                val[ptr[di] + i] = -1.0/one_ring_.size();
            }
            idx[ptr[di] + one_ring_.size()] = dim_*point_idx_+di;
            val[ptr[di] + one_ring_.size()] = 1.0;
        }
        return 0;
    }
    virtual size_t jac_nnz(void) const
    {
        return (one_ring_.size()+1)*dim_of_f();
    }
private:
    const std::vector<size_t> one_ring_;
    const size_t point_idx_;
    const size_t point_num_;
    const size_t dim_;
};

class fix_var_func : public jtf_func
{
public :
    fix_var_func(const size_t node_number,
                 const size_t idx0,
                 const double v)
        :node_number_(node_number), idx0_(idx0), v_(v) {}
    virtual ~fix_var_func() {}
    virtual size_t dim() const
    {
        return node_number_;
    }
    virtual int val(const val_type *x, val_type &v) const
    {
        v += 0.5*(x[idx0_] - v_)*(x[idx0_] - v_);
        return 0;
    }
    virtual int gra(const val_type *x, size_t &nnz, val_type *g, int_type *idx)
    {
        if(g == 0 && idx == 0) {
            nnz = 1;
            return 0;
        }

        g[0] = x[idx0_]-v_;
        idx[0] = idx0_;


        return 0;
    }
    virtual int gra(const val_type *x, val_type *g)
    {
        g[idx0_] += x[idx0_]-v_;

        return 0;
    }

    virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                    val_type *h, int_type *ptr, int_type *idx, val_type alpha)
    {
        if(h == 0  && idx == 0 && ptr == 0) {
            nnz = 1;
            format = 1;
            return 0;
        }
        if(h == 0 && idx != 0 && ptr != 0) {

            ptr[idx0_+1] = ptr[idx0_]+1;
            idx[ptr[idx0_]+0] = idx0_;


            return 0;
        }
        if(h != 0 && idx != 0 && ptr != 0) {
            jtf::function::add_to_csc(h,ptr,idx, idx0_, idx0_,1.0);

        }
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

template <typename val_type, typename int_type>
hj_func * build_degree_smooth_func(
    const zjucad::matrix::matrix<size_t> & mesh)
{
    using namespace std;
    using namespace zjucad::matrix;

    std::map<size_t,set<size_t> > point2neighbour;
    for(size_t ci = 0; ci < mesh.size(2); ++ci) {
        for(size_t pi = 0; pi < mesh.size(1); ++pi) {
            for(size_t pj = 1; pj < mesh.size(1); ++pj)
                point2neighbour[mesh(pi,ci)].insert(mesh((pj+pi)%mesh.size(1),ci));
        }
    }

    std::shared_ptr<vector<std::shared_ptr<const hj::function::function_t<val_type,int_type> > > >
    funcs(new vector<std::shared_ptr<const hj::function::function_t<val_type,int_type> > >);

    vector<size_t> neighbour;
    for(const auto & pi2n : point2neighbour) {
        neighbour.resize(pi2n.second.size());
        std::copy(pi2n.second.begin(), pi2n.second.end(), neighbour.begin());
        funcs->push_back(hj_func_cons_ptr(new degree_smooth_func(
                                              point2neighbour.size(), pi2n.first,neighbour)));
    }

    return hj::function::new_catenated_function<double, int_type>(funcs);
}

template <typename OS, typename FLOAT, typename INT>
void tet2vtk(
    OS &os,
    const FLOAT *node, size_t node_num,
    const INT *tet, size_t tet_num)
{
    os << "# vtk DataFile Version 2.0\nTET\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";
    os << "POINTS " << node_num << " float\n";
    for(size_t i = 0; i < node_num; ++i)
        os << node[i*3+0] << " " << node[i*3+1] << " " << node[i*3+2] << "\n";

    os << "CELLS " << tet_num << " " << tet_num*5 << "\n";
    for(size_t i = 0; i < tet_num; ++i)
        os << 4 << "  "
           << tet[i*4+0] << " " << tet[i*4+1] << " "
           << tet[i*4+2] << " " << tet[i*4+3] << "\n";
    os << "CELL_TYPES " << tet_num << "\n";
    for(size_t i = 0; i < tet_num; ++i)
        os << 10 << "\n";
}

void test_smoothing_soft_constraint(const char * tet_file)
{
    zjucad::matrix::matrix<size_t> tetmesh;
    zjucad::matrix::matrix<double> tet_node;
    if(jtf::mesh::tet_mesh_read_from_zjumat(tet_file, &tet_node, &tetmesh)) {
        std::cerr << "# can not open tet file." << std::endl;
        return;
    }

    zjucad::matrix::matrix<size_t> outside_surface;
    std::auto_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tetmesh));
    if(!fa.get()) {
        std::cerr << "# [error] can not build face2tet_adjacent." << std::endl;
        return;
    }
    jtf::mesh::get_outside_face(*fa, outside_surface);
    std::set<size_t> outside_points(outside_surface.begin(), outside_surface.end());

    hj_func_cons_ptr obj(build_degree_smooth_func<double, int32_t>(tetmesh));

    jtf_func_cons_ptr inner_obj(jtf::function::least_square_warpper(*obj));
    std::vector<jtf_func_cons_ptr>  all_func;
    {
        for(const std::set<size_t>::const_iterator cit = outside_points.begin();
                cit != outside_points.end(); ++cit) {
            all_func.push_back(
                jtf_func_cons_ptr(new node_fix2(tet_node.size(2), *cit,
                                                tet_node(zjucad::matrix::colon(),*cit),10)));
        }
        all_func.push_back(inner_obj);
    }

    jtf_func_cons_ptr total_obj(new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD>(all_func));

    ////////////////////////////////////////////////////////////////////////////////////////////////

    typedef kernel_traits<typename hj::sparse::csc<double,int32_t> > kernel_type;

    std::shared_ptr<func_container > eqn_func(new func_container),
        ineqn_func(new func_container);

    boost::shared_ptr<objective_function<kernel_type> > jftof(
        new jtf_func_to_obj_func<kernel_type>(total_obj, eqn_func, ineqn_func));


//    third_solver_ipopt tsb(jftof);
//    boost::property_tree::ptree pt;
//    pt.put("iter.value",100);
//    tsb.solve(tet_node, pt);
//    std::ofstream ofs("output.vtk");

//    tet2vtk(ofs, &tet_node[0], tet_node.size(2), &tetmesh[0], tetmesh.size(2));
}

void test_smoothing_hard_constraint(const char * tet_file)
{
    zjucad::matrix::matrix<size_t> tetmesh;
    zjucad::matrix::matrix<double> tet_node;
    if(jtf::mesh::tet_mesh_read_from_zjumat(tet_file, &tet_node, &tetmesh)) {
        std::cerr << "# can not open tet file." << std::endl;
        return;
    }

    zjucad::matrix::matrix<size_t> outside_surface;
    std::auto_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tetmesh));
    if(!fa.get()) {
        std::cerr << "# [error] can not build face2tet_adjacent." << std::endl;
        return;
    }
    jtf::mesh::get_outside_face(*fa, outside_surface);
    std::set<size_t> outside_points(outside_surface.begin(), outside_surface.end());

    hj_func_cons_ptr obj(build_degree_smooth_func<double, int32_t>(tetmesh));

    jtf_func_cons_ptr total_obj(jtf::function::least_square_warpper(*obj));

    ////////////////////////////////////////////////////////////////////////////////////////////////

    typedef kernel_traits<typename hj::sparse::csc<double,int32_t> > kernel_type;

    std::shared_ptr<func_container > eqn_func(new func_container),
        ineqn_func(new func_container);

    for(const std::set<size_t>::const_iterator cit = outside_points.begin();
            cit != outside_points.end(); ++cit) {
        for(size_t di =0; di < 3; ++di) {
            eqn_func->push_back(
                jtf_func_cons_ptr(new fix_var_func(tet_node.size(), 3*(*cit)+di,
                                                   tet_node(di,*cit))));
        }
    }

    boost::shared_ptr<objective_function<kernel_type> > jftof(
        new jtf_func_to_obj_func<kernel_type>(total_obj, eqn_func, ineqn_func));


//    third_solver_ipopt tsb(jftof);
//    boost::property_tree::ptree pt;
//    pt.put("iter.value",500);
//    tsb.solve(tet_node, pt);
//    std::ofstream ofs("output.vtk");

//    tet2vtk(ofs, &tet_node[0], tet_node.size(2), &tetmesh[0], tetmesh.size(2));
}
}
}
#endif // PROBLEM_SMOOTH_H
