#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <hjlib/sparse/sparse.h>
#include <hjlib/sparse/operation.h>
#include <fstream>
#include <functional>
#include <optimization.h>

#include "KKT.h"
#include "gmres.h"
#include "qmr.h"

#include "../../util/include/util.h"
#include <jtflib/function/operation.h>
#include <jtflib/function/func_aux.h>

#include "SQP.h"
using namespace std;
using namespace hj::sparse;
using namespace zjucad::matrix;


namespace jtf{
  jtf::function::functionN1_t<double,int32_t> *f4alglib = nullptr; // only used for alglib
  const KKT * kkt_ptr = nullptr; // only used for alglib

  int KKT::solve(double *x, boost::property_tree::ptree &pt)
  {
    pt.put("kkt_solver.desc", "bleic/jtf");
    const string kkt_solver = pt.get<string>("kkt_solver.value");
    if(kkt_solver == "bleic")
      return solve_alglib(x,pt);
    else if(kkt_solver == "jtf")
      return solve_sparse(x, pt);
  }

  /// ( H(x) A^T) *  (\delta x)      = -1*  (grad(x))
  /// ( A    0 )     (\delta \lamda)        (Ax-b)
  int KKT::solve_sparse(double *x, boost::property_tree::ptree &pt)
  {
    //preprocess_init_data(x,pt);

    const double epsf = pt.get<double>("epsf.value",1e-6);
    const double epsg = pt.get<double>("epsg.value",1e-6);
    const double epsx = pt.get<double>("epsx.value",1e-6);
    const double epsl = pt.get<double>("epsl.value",1e-6);
    const size_t max_iter = pt.get<size_t>("iter.value",1000);

    assert(f_ && eqn_cons_);
    itr_matrix<double*> x_m(f_->dim(),1,x);

    hj::sparse::csc<double,int32_t> H;
    size_t hnnz = 0,format;
    {// get Hessian of f
      f_->hes(x, hnnz, format, 0,0,0);
      H.resize(f_->dim(),f_->dim(), hnnz);
      f_->hes(x, hnnz, format, 0, &H.ptr()[0],&H.idx()[0]);
    }

    hj::sparse::csc<double,int32_t> A,AT;
    matrix<double> b;
    get_A_b(A,b);
    trans(A,AT);

    zjucad::matrix::matrix<double> grad(f_->dim(),1);

    matrix<double> Ax_b(b.size(1),1);
    matrix<double> grad_Ax_b(f_->dim() + b.size(1),1);
    //matrix<double> K = zeros<double>(f_->dim() + b.size(1), f_->dim() + b.size(1));
    hj::sparse::csc<double,int32_t> K(f_->dim() + b.size(1),
                                      f_->dim() + b.size(1),
                                      hnnz + 2 * hj::sparse::nnz(A));

    {// assemble K
      assert(A.ptr().size() == H.ptr().size());
      // assemble left two blocks
      for(size_t pi = 0; pi < A.ptr().size()-1; ++pi){
          K.ptr()[pi+1] = K.ptr()[pi] + H.ptr()[pi+1] - H.ptr()[pi]
              + A.ptr()[pi+1] - A.ptr()[pi];
          for(size_t i = H.ptr()[pi]; i < H.ptr()[pi+1]; ++i){
              K.idx()[K.ptr()[pi] + i - H.ptr()[pi]] = H.idx()[i];
              K.val()[K.ptr()[pi] + i - H.ptr()[pi]] = H.val()[i];
            }
          for(size_t i = A.ptr()[pi]; i < A.ptr()[pi+1]; ++i){
              int32_t idx = K.ptr()[pi] + i - A.ptr()[pi] + H.ptr()[pi+1]- H.ptr()[pi];
              K.idx()[idx]
                  = A.idx()[i] + A.size(2);
              K.val()[idx] = A.val()[i];
            }
        }
      for(size_t pi = 0; pi < AT.ptr().size()-1; ++pi){
          K.ptr()[pi+1 + A.size(2)] = K.ptr()[pi + A.size(2)] + AT.ptr()[pi+1] - AT.ptr()[pi];
          for(size_t i = AT.ptr()[pi]; i < AT.ptr()[pi+1]; ++i){
              K.idx()[K.ptr()[pi +  A.size(2)] + i - AT.ptr()[pi] ] = AT.idx()[i];
              K.val()[K.ptr()[pi +  A.size(2)] + i - AT.ptr()[pi] ] = AT.val()[i];
            }
        }
    }

    matrix<double> dx_dlambda = zeros<double>(f_->dim() + b.size(1),1);
    itr_matrix<double*> dx(f_->dim(), 1, &dx_dlambda[0]);
    itr_matrix<double*> dlambda(b.size(1),1, &dx_dlambda[f_->dim()]);
    matrix<double> x_temp = x_m;
    matrix<double> lambda = ones<double>(b.size(),1);
    matrix<double> lambda_temp = lambda;

    matrix<double> r_temp = zeros<double>(x_temp.size() + lambda_temp.size(),1);
    matrix<double> Ax_b_temp = zeros<double>(b.size(),1);

    double f;
    double gn = 0;

    static jtf::PETsc_init pi;
    for(size_t i = 0; i < max_iter; ++i){
        grad *= 0;
        H.val() *= 0;
        f_->hes(&x[0], hnnz, format, &H.val()[0], &H.ptr()[0], &H.idx()[0]);
        f_->gra(&x[0],&grad[0]);

        Ax_b = -1*b;

        hj::sparse::mv(false, A, x_m, Ax_b); // Ax_b += A*x
        hj::sparse::mv(false, AT, lambda , grad); // grad += AT*lambda

        std::copy(grad.begin(), grad.end(), grad_Ax_b.begin());
        std::copy(Ax_b.begin(), Ax_b.end(), grad_Ax_b.begin() + grad.size());
        grad_Ax_b *= -1;

        {// assign H to K block
          for(size_t pi = 0; pi < H.ptr().size()-1; ++pi){
              for(size_t i = H.ptr()[pi]; i < H.ptr()[pi+1]; ++i){
                  K.val()[K.ptr()[pi] + i - H.ptr()[pi]] = H.val()[i];
                }
            }
        }

        //shared_ptr<jacobi_precondition_sparse> jp(new jacobi_precondition_sparse(K,grad_Ax_b));
        //gmres_solver_sparse gs(K,grad_Ax_b);
        //qmr_solver_sparse qs(K, grad_Ax_b);

        PETsc_TFQMR qs(false, K, grad_Ax_b);
        //qs.solve2(&dx_dlambda[0], 1000, grad_Ax_b.size(), jp);

        //std::function<int(void)> func = std::bind(&qmr_solver_sparse::solve, &qs, &dx_dlambda[0], f_->dim());
        //cerr << "time " << jtf::util::func_time(func)  << "(s)." << endl;
        int rtn = qs.solve(&dx_dlambda[0], f_->dim());

        double t = 1.1;
#if 1
        {// backtracking
          double alpha = 5e-3;
          double beta = 0.9;
          const double norm_grad_Ax_b = norm(grad_Ax_b);

          while(1){
              x_temp = x_m + t * dx;
              lambda_temp = lambda + t * dlambda;
              Ax_b_temp = -1.0*b;
              grad *= 0;
              f_->gra(&x_temp[0], &grad[0]);
              hj::sparse::mv(false, AT, lambda_temp, grad);
              hj::sparse::mv(false, A, x_temp, Ax_b_temp);

              std::copy(grad.begin(), grad.end(), r_temp.begin());
              std::copy(Ax_b_temp.begin(), Ax_b_temp.end(), r_temp.begin() + x_m.size());

              double new_norm = norm(r_temp);
              if(new_norm > (1-alpha * t) * norm_grad_Ax_b) t *= beta;
              else {
                  if(fabs(t) < 1e-7) {// t ==0
                      t = 1.1;
                      alpha *= 0.618;
                      if(fabs(alpha) < 1e-7) {
                          t = 0;
                          break;
                        }
                    }else break;
                };
            }
        }
#else
        {
          double t1 = 0, t2 = 100;
          double t = 1;
          double rho = 0.3;

          matrix<double> Ax_b;
          double varphi_0 = 0;
          {
            f_->val(x, varphi_0);
            Ax_b = -1 * b;
            hj::sparse::mv(false, A, x_m, Ax_b);
            varphi_0 += dot(lambda,  Ax_b);
          }
          matrix<double> partial_varphi_0 = grad_Ax_b;


          while(1){
              x_temp = x_m + t * dx;
              lambda_temp = lambda + t * dlambda;

              double varphi_t = 0;
              {
                f_->val(&x_temp[0], varphi_t);
                Ax_b = -1 * b;
                hj::sparse::mv(false, A, x_temp, Ax_b);
                varphi_t += dot(lambda_temp, Ax_b);
              }

              itr_matrix<double*> pv_x(f_->dim(),1, &partial_varphi_0[0]);
              itr_matrix<double*> pv_lambda(dlambda.size(),1, &partial_varphi_0[0]+f_->dim());
              const double dot_tdk_partial = t*(dot(dx, pv_x) + dot(dlambda,pv_lambda));
              if(varphi_t > varphi_0 + rho * dot_tdk_partial){
                  t2 = t;
                  t = (t1+t2)/2;
                  continue;
                }else{
                  if(varphi_t < varphi_0 + (1-rho)*dot_tdk_partial){
                      t1 = t;
                      if(fabs(t2-100) < 1e-7) t = 0.9*t;
                      else t = (t1+t2)/2;
                      continue;
                    }else break;
                }
            }
        }
#endif
        dx_dlambda *= t;
        x_m += dx; // need a linear search
        lambda += dlambda;
        f = 0;
        f_->val(x, f);

        Ax_b_temp = b * -1.0;
        hj::sparse::mv(false, A,x_m, Ax_b_temp);
        gn = norm(grad);
        cerr << "# " << i << "\t f " << f << "\t g " << gn  << "\t Ax-b " << norm(Ax_b_temp)<< endl;
        if((fabs(f) < epsf || gn < epsg || norm(x_m) < epsx ) &&
           norm(Ax_b_temp) < epsl)
          return 0;
      }

    return __LINE__;
  }


  void function_grad(const alglib::real_1d_array &x, double &func,
                     alglib::real_1d_array &grad, void *ptr)
  {
    assert(f4alglib != nullptr);
    assert(kkt_ptr);

    static hj::sparse::csc<double,int32_t> A;
    static zjucad::matrix::matrix<double> b;

    if(A.size(1) == 0 || b.size() == 0){
        kkt_ptr->get_A_b(A,b);
      }

    const double *x_raw = x.getcontent();
    double *grad_raw = grad.getcontent();
    itr_matrix<const double*> x_m(f4alglib->dim(),1,x_raw);
    itr_matrix<double*> grad_m(f4alglib->dim(),1,grad_raw);

    func = 0;
    grad_m *= 0;
    f4alglib->val(x_raw, func);
    f4alglib->gra(x_raw, grad_raw);

    matrix<double> r = -1 * b;
    hj::sparse::mv(false, A, x_m, r);

    cerr << "# [info] \t f " << func << "\t g " << norm(grad_m)
         << "\t Ax-b " << norm(r) << endl;
  }


  int KKT::solve_alglib(double *x, boost::property_tree::ptree &pt)
  {
    using namespace alglib;

    minbleicstate state;
    minbleicreport rep;

    const double epsf = pt.get<double>("epsf.value",1e-6);
    const double epsg = pt.get<double>("epsg.value",1e-6);
    const double epsx = pt.get<double>("epsx.value",1e-6);
    const double epsl = pt.get<double>("epsl.value",1e-6);
    ae_int_t max_its = pt.get<ae_int_t>("iter.value",1000);

    real_1d_array ra;
    ra.setcontent(f_->dim(), x);

    real_2d_array c;
    integer_1d_array ct; // determine constraint type 0 for =

    hj::sparse::csc<double,int32_t> A;
    matrix<double> b;
    get_A_b(A,b);

    {
      matrix<double> A_m;
      hj::sparse::convert(A,A_m);
      matrix<double> A_m_b(A_m.size(1), A_m.size(2)+1);
      A_m_b(colon(),colon(0,A_m.size(2)-1)) = A_m;
      A_m_b(colon(),A_m.size(2)) = b;

      c.setcontent(A_m_b.size(1), A_m_b.size(2), &A_m_b[0]);
    }


    zjucad::matrix::matrix<ae_int_t> c_m = zeros<ae_int_t>(b.size(),1);
    ct.setcontent(c_m.size(), &c_m[0]);

    minbleiccreate(ra, state);
    minbleicsetlc(state, c, ct);
    //    minbleicsetcond(state, epsg, epsf, epsx, max_its);

    f4alglib = f_;
    kkt_ptr = this;
    minbleicoptimize(state, function_grad);
    minbleicresults(state, ra, rep);

    double *x_finial = ra.getcontent();
    std::copy(x_finial, x_finial+f_->dim(), x);
    int rtn = rep.terminationtype;
    switch (rtn) {
      case -7 : {
          cerr <<  "gradient verification failed." << endl; break;
        }
      case -3 : {
          cerr <<  "inconsistent constraints. Feasible point is"
                   "either nonexistent or too hard to find. Try to"
                   "restart optimizer with better initial approximation" << endl;
          break;
        }
      case 1 : {
          cerr << "relative function improvement is no more than EpsF." << endl;
          break;
        }
      case 2 : {cerr << "scaled step is no more than EpsX." << endl; break;}
      case 4 : { cerr << "scaled gradient norm is no more than EpsG." << endl; break;}
      case 5 : { cerr << "MaxIts steps was taken." << endl; break;}

      default:
        break;
      }

    return 0;
  }
  int KKT::get_A_b(hj::sparse::csc<double,int32_t> & A,
                   matrix<double> & b) const
  {
    assert(eqn_cons_);
    matrix<double> x_temp = zeros<double>(f_->dim(),1);
    b = zeros<double>(eqn_cons_->size(),1);
    size_t nnz = 0;
    size_t total_nnz = 0;
    for(size_t fi = 0; fi < eqn_cons_->size(); ++fi){
        const_cast<jtf::function::functionN1_t<double,int32_t> *>((*eqn_cons_)[fi].get())->gra(&x_temp[0],nnz, 0,0);
        //			x_temp.resize(1,f_->dim(), nnz, &x_temp.val()[0], &x_temp.ptr()[0], &x_temp.idx()[0]);
        total_nnz += nnz;
        (*eqn_cons_)[fi]->val(&x_temp[0], b[fi]);
      }
    b *= -1;

    // it's AT;
    hj::sparse::csc<double,int32_t> AT(f_->dim(), eqn_cons_->size(),total_nnz);

    matrix<double> g(total_nnz,1);
    matrix<int32_t> idx(total_nnz, 1);

    for(size_t fi = 0; fi < eqn_cons_->size(); ++fi){
        const_cast<jtf::function::functionN1_t<double,int32_t> *>((*eqn_cons_)[fi].get())->gra(&x_temp[0], nnz, 0, 0);
        const_cast<jtf::function::functionN1_t<double,int32_t> *>((*eqn_cons_)[fi].get())->gra(&x_temp[0], nnz, &g[0], &idx[0]);
        AT.ptr()[fi+1] = AT.ptr()[fi] + nnz;
        std::copy(g.begin(), g.begin() + nnz, &AT.val()[AT.ptr()[fi]]);
        std::copy(idx.begin(), idx.begin() + nnz, &AT.idx()[AT.ptr()[fi]]);
      }

    trans(AT, A);

    return 0;
  }

//  int KKT::preprocess_init_data(double *x,boost::property_tree::ptree &pt)
//  {

//    if(ineqn_cons_ == nullptr) return 0;

//    vector<jtf::function::functionN1_t<double,int32_t> const* > obj;
//    obj.resize(ineqn_cons_->size() + 1);

//    obj[0] = f_;

//    vector<shared_ptr<const jtf::function::functionN1_t<double,int32_t> > > ineqn_obj;
//    double v = 0;
//    while(1){
//        ineqn_obj.clear();
//        size_t vio_ineqn_num = 0;
//        //for(size_t fi = 0; fi < ineqn_cons_->size(); ++fi){
//        for(size_t fi = 0; fi < ineqn_cons_->size(); ++fi){
//            v = 0;
//            (*ineqn_cons_)[fi]->val(x,v);
//            if(v < 0){
//                ++vio_ineqn_num;
//                ineqn_obj.push_back(
//                      std::shared_ptr<const jtf::function::functionN1_t<double,int32_t> >(
//                        (*ineqn_cons_)[fi]*-1000.0));
//                obj[fi+1] =ineqn_obj.back().get();
//              }else{
//                ineqn_obj.push_back(std::shared_ptr<const jtf::function::functionN1_t<double,int32_t> >(
//                                      jtf::function::pow_warpper((*ineqn_cons_)[fi],-1.0)));


//                obj[fi+1] = ineqn_obj.back().get();
//              }
//          }

//        cerr << "# [info] prprocess init x, violate ineqn: " << vio_ineqn_num
//             << " satisfy ineqn: " << ineqn_cons_->size() - vio_ineqn_num << endl;
//        if(vio_ineqn_num == 0) return 0;

//        shared_ptr<jtf::function::functionN1_t<double,int32_t> > all_func(
//              new jtf::function::sum_function<double,int32_t, jtf::function::RAW>(obj));

//        jtf::SQP sqp;
//        sqp.set_f(*all_func);

//        sqp.solve(x, pt);
//      }
//    return __LINE__;
//  }

}
