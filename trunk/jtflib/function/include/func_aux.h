#ifndef JTF_FUNCTION_AUX_H_
#define JTF_FUNCTION_AUX_H_

#include <vector>
#include <iostream>

#include <hjlib/function/function.h>
#include <zjucad/matrix/matrix.h>

#include "function.h"

namespace jtf { namespace function {

    template <typename val_type, typename int_type0, typename int_type1>
    int add_to_csc(val_type *val, const int_type0 *ptr, const int_type0 *idx,
                   const int_type1 r, const int_type1 c, val_type v) {
      size_t off;
      for(off = ptr[c]; off < ptr[c+1]; ++off) {
          if(idx[off] == r)
            break;
        }
      if(off == ptr[c+1]) {
          std::cerr << "# incorrect sparse pattern in add_to_csc: "
                    << r << ":" << c;
          for(off = ptr[c]; off < ptr[c+1]; ++off) {
              std::cerr << " " << idx[off];
            }
          std::cerr << std::endl;
          return __LINE__;
        }
      val[off] += v;
      return 0;
    }


    template <typename val_type, typename int_type>
    double gra_err(jtf::function::functionN1_t<val_type,int_type> &f, val_type *x)
    {
      using namespace zjucad::matrix;
      const size_t dim = f.dim();
      double val = 0;
      f.val(x, val);
      matrix<val_type> g = zjucad::matrix::zeros<val_type>(dim, 1);
      f.gra(x, &g[0]);
      std::cerr << "max g: " << zjucad::matrix::max(zjucad::matrix::fabs(g)) << std::endl;
      const val_type eps = 1e-6;
      for(size_t xi = 0; xi < dim; ++xi) {
          const val_type save = x[xi];
          val_type v[2] = {0, 0};
          x[xi] = save-eps;
          f.val(x, v[0]);
          x[xi] = save+eps;
          f.val(x, v[1]);
          g[xi] -= (v[1]-v[0])/(2*eps);
          x[xi] = save;
        }
      return zjucad::matrix::max(zjucad::matrix::fabs(g));
    }

    template <typename val_type, typename int_type>
    double hes_err(jtf::function::functionN1_t<val_type,int_type> &f, val_type *x)
    {
      using namespace zjucad::matrix;
      const size_t dim = f.dim();
      matrix<val_type> hes(dim, dim);
      size_t nnz, format = -1;
      f.hes(x, nnz, format, 0, 0, 0);
      matrix<val_type> h = zeros<val_type>(nnz, 1);
      matrix<int_type> ptr(dim+1), idx(nnz);
      //TODO: only accept csc format
      assert(format == 1);
      f.hes(x, nnz, format, 0, &ptr[0], &idx[0]);
      f.hes(x, nnz, format, &h[0], &ptr[0], &idx[0]);
      for(size_t ci = 0; ci < dim; ++ci) {
          for(size_t nzi = ptr[ci]; nzi < ptr[ci+1]; ++nzi)
            hes(idx[nzi], ci) = h[nzi];
        }
      std::cout << "max hes: " << zjucad::matrix::max(zjucad::matrix::fabs(hes)) << std::endl;

      const val_type eps = 1e-6;
      matrix<val_type> g0 = zjucad::matrix::zeros<val_type>(dim, 1),
          ga = zjucad::matrix::zeros<val_type>(dim, 1),
          gb = zjucad::matrix::zeros<val_type>(dim, 1);

      f.gra(x, &g0[0]);
      for(size_t xi = 0; xi < dim; ++xi) {
          const val_type x0 = x[xi];

          x[xi] = x0+eps;
          ga(zjucad::matrix::colon()) = 0;
          f.gra(x, &ga[0]);

          x[xi] = x0-eps;
          gb(zjucad::matrix::colon()) = 0;
          f.gra(x, &gb[0]);

          hes(zjucad::matrix::colon(), xi) -= (ga-gb)/(2*eps);

          x[xi] = x0;
        }
      return zjucad::matrix::max(zjucad::matrix::fabs(hes));
    }

    ///
    /// @brief least_square_warpper
    /// @param f input instance of function_t, you should ensure f is alive during the whole processing
    /// @return
    ///
    jtf::function::functionN1_t<double,int32_t> *
    least_square_warpper(const hj::function::function_t<double,int32_t> & f);

    jtf::function::functionN1_t<double,int32_t> *
    least_square_warpper(const std::shared_ptr<const hj::function::function_t<double,int32_t> > fptr);

    ///
    /// \brief neg_log_warpper
    /// \param f
    /// \return
    ///
    jtf::function::functionN1_t<double,int32_t> *
    neg_log_warpper(jtf::function::functionN1_t<double,int32_t> &f);

    jtf::function::functionN1_t<double,int32_t> *
    neg_log_warpper(std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > fptr);

    jtf::function::functionN1_t<double,int32_t> *
    neg_log_warpper(std::shared_ptr<hj::function::function_t<double,int32_t> > hptr,
                    const zjucad::matrix::matrix<double> & weight);


    ///
    /// \brief function convertions
    /// \param func
    /// \return
    ///
    jtf::function::functionN1_t<double,int32_t> *
    hj_func_to_jtf_func(const hj::function::function_t<double,int32_t> * func);

    hj::function::function_t<double,int32_t> *
    jtf_func_to_hj_func(const jtf::function::functionN1_t<double,int32_t> *func);


//    ///
//    /// \brief pow_warpper
//    /// \param func
//    /// \return
//    ///
//    jtf::function::functionN1_t<double,int32_t> *
//    pow_warpper(const jtf::function::functionN1_t<double,int32_t>&func,double n);

//    jtf::function::functionN1_t<double,int32_t> *
//    pow_warpper(const std::shared_ptr<const jtf::function::functionN1_t<double,int32_t> > func,
//                double n);

  }}

#endif
