#ifndef PROBLEM_FRAME_ALIGN_H
#define PROBLEM_FRAME_ALIGN_H


#include <jtflib/function/function.h>
#include <jtflib/function/func_aux.h>
#include <hjlib/sparse/sparse.h>
#include <hjlib/function/function.h>
#include <zjucad/matrix/itr_matrix.h>

#include "jtf_func_to_obj_func.h"
#include "third_solver_ipopt.h"

namespace  zjucad
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

template<typename val_type, typename int_type>
class frame_func : public hj::function::function_t<val_type,int_type>
{
public:
    frame_func(const zjucad::matrix::matrix<val_type> & node,
               const zjucad::matrix::matrix<size_t> & tet,
               const size_t idx,
               const zjucad::matrix::matrix<val_type> &M_inv,
               const zjucad::matrix::matrix<double> &frame,
               const val_type w)
        :grad_op_(M_inv), tet_(tet), node_num_(node.size(2)), idx_(idx) , weight_(w),
         nnz_(-1), frame_(trans(frame))
    {
        using namespace zjucad::matrix;
        one_tet_ = tet_(colon(),idx_);
    }
    virtual ~frame_func() {}
public:

    virtual size_t dim_of_x(void) const
    {
        return node_num_*3;
    }
    virtual size_t dim_of_f(void) const
    {
        return 9;
    }
    virtual int val(const val_type *x, val_type *f,
                    hj::function::func_ctx *ctx = 0) const
    {
        using namespace zjucad::matrix;
        zjucad::matrix::itr_matrix<val_type *> f0(3, 3, f);

        const zjucad::matrix::itr_matrix<const val_type *> T(3, node_num_, x);
        f0 = T(colon(), one_tet_)*grad_op_;

        f0 -= frame_;
        f0 *= weight_;

        //for(size_t i = 0; i < f0.size(); ++i) erase_nan(f0[i]);
        return 0;
    }
    virtual int jac(const val_type *x, val_type *val, int_type *ptr = 0,
                    int_type *idx = 0, hj::function::func_ctx *ctx = 0) const
    {
        for(int c = 0; c < 3; ++c) {
            for(int r = 0; r < 3; ++r) {
                const int fi = c*3+r;
                ptr[fi+1] = ptr[fi] + tet_.size(1);
                for(int k = 0; k < tet_.size(1); ++k) {
                    idx[ptr[fi]+k] = tet_(k, idx_)*3+r;
                    val[ptr[fi]+k] = grad_op_(k, c) * weight_;
                    // erase_nan(val[ptr[fi]+k]);
                }
            }
        }

        return 0;
    }
    virtual size_t jac_nnz(void) const
    {
        return 4*dim_of_f();
    }

private:
    const size_t idx_, node_num_;
    const zjucad::matrix::matrix<size_t> &tet_;
    zjucad::matrix::matrix<val_type> grad_op_;
    const double weight_;
    zjucad::matrix::matrix<size_t> one_tet_;
    const zjucad::matrix::matrix<double> frame_;
    size_t nnz_ ;
};

template <typename E>
inline typename E::value_type cal_tet_vol(const zjucad::matrix::matrix_expression<E> &tet)
{
    using namespace zjucad::matrix;
    assert(tet().size(1) == 3);
    if(tet().size(2) == 4) {
        zjucad::matrix::matrix<typename E::value_type> edges
            = tet()(colon(), colon(1, 3)) - tet()(colon(), 0)*
              ones<typename E::value_type>(1, 3);
        return cal_tet_vol(edges);
    }
    if(tet().size(2) == 3) // fix by jtf, here tet means three edges (3*3)
        return dot(cross(tet()(colon(), 0), tet()(colon(), 1)), tet()(colon(), 2))/6.0;
}

template <typename val_type, typename int_type>
hj_func * build_frame_func(
    const zjucad::matrix::matrix<size_t> & mesh,
    const zjucad::matrix::matrix<val_type> & node,
    const zjucad::matrix::matrix<zjucad::matrix::matrix<val_type> > &M_inv,
    const zjucad::matrix::matrix<zjucad::matrix::matrix<val_type> > &frame)
{
    using namespace std;
    using namespace zjucad::matrix;

    zjucad::matrix::matrix<val_type> weight = ones<val_type>(mesh.size(2),1);

    zjucad::matrix::matrix<val_type> one_simplex_node = zeros<val_type>(3,mesh.size(1));
    for(size_t ti = 0; ti < mesh.size(2); ++ti) {
        one_simplex_node = node(colon(), mesh(colon(),ti));
        weight[ti] = std::fabs(cal_tet_vol(one_simplex_node));
    }

    const double total_w = std::accumulate(weight.begin(), weight.end(), 0.0);

    std::shared_ptr<vector<std::shared_ptr<hj::function::function_t<val_type,int_type> > > >
    funcs(new vector<std::shared_ptr<hj::function::function_t<val_type,int_type> > >(mesh.size(2)));

    for(size_t ti = 0; ti < mesh.size(2); ++ti) {
        (*funcs)[ti].reset(
            new frame_func<val_type, int_type>(
                node, mesh, ti, M_inv[ti], frame[ti], std::sqrt(weight[ti]/total_w)));
    }

    return hj::function::new_catenated_function<val_type, int_type>(funcs);
}

class fix_var_func4 : public jtf_func
{
public :
    fix_var_func4(const size_t node_number,
                  const size_t idx0,
                  const double v)
        :node_number_(node_number), idx0_(idx0), v_(v) {}
    virtual ~fix_var_func4() {}
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

class fix_var_func5 : public jtf_func
{
public :
    fix_var_func5(const size_t node_number,
                  const size_t idx0,
                  const double v)
        :node_number_(node_number), idx0_(idx0), v_(v) {}
    virtual ~fix_var_func5() {}
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

template <typename OS, typename FLOAT, typename INT>
void tet2vtk4(
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

int tet_mesh_read_from_vtk2( const char *path,
                             zjucad::matrix::matrix<double> &node,
                             zjucad::matrix::matrix<size_t> &tet )
{
    std::ifstream ifs(path);
    if(ifs.fail()) {
        std::cerr << "[info] " << "can not open file" << path << std::endl;
        return __LINE__;
    }

    std::string str;
    int point_num = 0,cell_num = 0;

    while(!ifs.eof()) {
        ifs >> str;
        if(str == "POINTS") {
            ifs >> point_num >> str;
            node.resize(3, point_num);
            for(size_t i = 0; i < point_num; ++i) {
                for(size_t j = 0; j < 3; ++j)
                    ifs >> node(j, i);
            }
            continue;
        }
        if(str == "CELLS") {
            ifs >> cell_num >> str;
            int point_number_of_cell = 0;
            std::vector<size_t> tet_temp;
            for(size_t ci = 0; ci < cell_num; ++ci) {
                ifs >> point_number_of_cell;
                if(point_number_of_cell != 4) {
                    for(size_t i = 0; i < point_number_of_cell; ++i)
                        ifs >> str;
                } else {
                    int p;
                    for(size_t i = 0; i < point_number_of_cell; ++i) {
                        ifs >> p;
                        tet_temp.push_back(p);
                    }
                }
            }
            tet.resize(4, tet_temp.size()/4);
            std::copy(tet_temp.begin(), tet_temp.end(), tet.begin());
        }
    }
    return 0;
}

int read_M_inv2(const char * file,
                zjucad::matrix::matrix<zjucad::matrix::matrix<double> > &M_inv)
{
    std::ifstream ifs(file);
    if(ifs.fail()) {
        std::cerr << "can not load read_M_inv" << std::endl;
        return __LINE__;
    }
    size_t num;
    ifs >> num;
    M_inv.resize(num,1);
    zjucad::matrix::matrix<double> one_inv(4,3);
    for(size_t i = 0 ; i < num; ++i) {
        for(size_t j = 0; j < 12; ++j) {
            ifs >> one_inv[j];
        }
        M_inv[i] = one_inv;
    }
    return 0;
}

int read_frame(const char *frame_file,
               zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & frame)
{
    std::ifstream ifs(frame_file);
    if(ifs.fail()) {
        std::cerr << "# [error] can not open frame file" << std::endl;
        return __LINE__;
    }
    size_t num;
    ifs >> num;
    frame.resize(num,1);
    for(size_t i = 0 ; i <num; ++i) {
        frame[i].resize(3,3);
        for(size_t j = 0; j < 9; ++j)
            ifs >> frame[i][j];
    }

    return 0;
}

void test_align_frame(const char * tet_file, const char * frame_file, const char * M_inv_file)
{
    zjucad::matrix::matrix<size_t> tetmesh;
    zjucad::matrix::matrix<double> tet_node;
    if(tet_mesh_read_from_vtk2(tet_file, tet_node, tetmesh))
        return __LINE__;

    zjucad::matrix::matrix<zjucad::matrix::matrix<double> > M_inv, frame;
    read_M_inv2(M_inv_file, M_inv);
    read_frame(frame_file, frame);

    hj_func_cons_ptr hptr(build_frame_func<double,int32_t>(tetmesh, tet_node,M_inv,frame));

    jtf_func_cons_ptr total_obj(jtf::function::least_square_warpper(*hptr));


    //////////////////////////////////////////////////////////////////////////////////////////////

    typedef kernel_traits<typename hj::sparse::csc<double,int32_t> > kernel_type;

    std::shared_ptr<func_container> eqn_func(new func_container),
        ineqn_func(new func_container);

    {
        for(size_t i = 0; i < tet_node.size(2); ++i) {
            if(tet_node(1,i) < 0) {
                eqn_func->push_back(jtf_func_cons_ptr(new fix_var_func4(tet_node.size(2), 3*i+1, 0.0)));
            } else {
                ineqn_func->push_back(jtf_func_cons_ptr(new fix_var_func5(tet_node.size(2), 3*i+1, 1.0)));
            }
        }
    }

    boost::shared_ptr<objective_function<kernel_type> >jftof(
        new jtf_func_to_obj_func<kernel_type>(total_obj,eqn_func, ineqn_func));

    third_solver_ipopt tsb(jftof);
    boost::property_tree::ptree pt;
    pt.put("iter.value",500);
    tsb.solve(tet_node, pt);
    {
        std::ofstream ofs("opt_tet.vtk");
        tet2vtk4(ofs, &tet_node[0], tet_node.size(2), &tetmesh[0], tetmesh.size(2));
    }
}
}
}

#endif // PROBLEM_FRAME_ALIGN_H
