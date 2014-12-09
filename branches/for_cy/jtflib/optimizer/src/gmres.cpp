#include "gmres.h"
#include <iostream>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>

using namespace std;
using namespace zjucad::matrix;

namespace jtf{
  typedef matrix<double> matrixd;

  template<class Real>
  void generate_plane_rotation(const Real &dx, const Real &dy, Real &cs, Real &sn)
  {
    if (fabs(dy) < 1e-12) {
        cs = 1.0;
        sn = 0.0;
      } else if (fabs(dy) > fabs(dx)) {
        Real temp = dx / dy;
        sn = 1.0 / sqrt( 1.0 + temp*temp );
        cs = temp * sn;
      } else {
        Real temp = dy / dx;
        cs = 1.0 / sqrt( 1.0 + temp*temp );
        sn = temp * cs;
      }
  }


  template<class Real>
  void apply_plane_rotation(Real &dx, Real &dy, const Real &cs, const Real &sn)
  {
    Real temp  =  cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = temp;
  }

  template<typename T>
  void update(zjucad::matrix::matrix_expression<T> &x, int k, const vector<matrixd> &H,
              const matrixd &s, const vector<matrixd> &v)
  {
    matrixd y = s;

    // Backsolve:
    for (int i = k; i >= 0; i--) {
        y[i] /= H[i][i];
        for (int j = i - 1; j >= 0; j--)
          y[j] -= H[i][j] * y[i];
      }

    for (int j = 0; j <= k; j++)
      x() += v[j] * y[j];
  }


  static matrixd operator *(const vector<matrixd> &A, const matrixd &b)
  {
    assert(A.size() == b.size(1));
    itr_matrix<const double*> A_mat(A[0].size(), A.size(), &A[0][0]);

    matrixd tmp = A_mat;
    tmp = temp(tmp * b);
    return tmp;
  }

  int gmres_solver::solve(double *x, size_t m, size_t max_iter) const
  {
    itr_matrix<double*> x_m(b_.size(),1,x);
    double resid = 0;
    size_t j = 0;
    matrixd s = zeros<double>(m+1,1);
    matrixd cs(m+1,1), sn(m+1,1);
    matrixd w,r = b_ - A_* x_m;
    double rho = norm(r);
    s[0] = rho;

    vector<matrixd> v(m+1);
    matrixd hk(m+1,1);
    vector<matrixd> H;
    H.reserve(m+1);
    size_t i;
    while(j < max_iter){
        v[0] = r/rho;
        s *= 0;
        s[0] = rho;

        for(i = 0; i < m && j <= max_iter; ++i, ++j){
            w = A_ * v[i];
            hk *=0; // clear hk
            for(size_t k = 0; k <= i;++k){
                hk[k] = dot(w, v[k]);
                w -= hk[k] * v[k];
              }
            hk[i+1] = norm(w);
            if(fabs(hk[i+1])<1e-7)
              v[i+1] = w;
            else
              v[i+1] = w / hk[i+1];  //TODO: how to address hk[i+1] == 0 case?

            {//Given rotations
              for(size_t k = 0; k < i; ++k)
                apply_plane_rotation(hk[k], hk[k+1], cs[k], sn[k]);

              generate_plane_rotation(hk[i],hk[i+1],cs[i],sn[i]);
              apply_plane_rotation(hk[i],hk[i+1],cs[i],sn[i]);
              apply_plane_rotation(s[i],s[i+1],cs[i],sn[i]);

              H.push_back(hk);
            }

            if((resid = fabs(s[i+1]/norm(b_))) < eps_){
                update(x_m,i,H,s,v);
                cerr << "# iter " << j << " " << resid << endl;
                cerr << "# converge." << endl;
                return 0;
              }
            cerr << "# iter " << j << " " << resid << endl;
          }

        update(x_m,i-1,H,s,v);
        r = b_ - A_ * x_m;
        rho = norm(r);
        if(fabs(rho/norm(b_)) < eps_){
            cerr << "# restart " << rho << endl;
            cerr << "# converge." << endl;
            return 0;
          }
        H.clear();
      }

    cerr << "# exceed max iterations." << endl;
    return 1;
  }


  int gmres_solver::solve2(double *x, size_t m, size_t max_iter,
                           const shared_ptr<jacobi_precondition> pre)const
  {
    itr_matrix<double*> x_m(b_.size(),1,x);
    double resid = 0;
    size_t j = 0;
    matrixd s = zeros<double>(m+1,1);
    matrixd cs(m+1,1), sn(m+1,1);
    matrixd w = b_,r = b_;

    matrixd mb = b_;
    pre->get_b(&mb[0]);
    double normb = norm(mb);
    pre->get_r(x, &r[0]); // M*(b-Ax)

    double rho = norm(r);
    s[0] = rho;

    if(fabs(normb) < 1e-7)
      normb = 1.0;

    if((resid = rho/normb) < eps_){
        cerr << "# converge." << endl;
        return 0;
      }

    vector<matrixd> v(m+1);
    matrixd hk(m+1,1);
    vector<matrixd> H(m+1);
    size_t i;
    while(j <= max_iter){
        v[0] = r/rho;
        s *= 0;
        s[0] = rho;

        for(i = 0; i < m && j <= max_iter; ++i, ++j){
            //w = A_ * v[i];
            pre->get_Ax(&v[i][0], &w[0]);
            hk *=0; // clear hk
            for(size_t k = 0; k <= i;++k){
                hk[k] = dot(w, v[k]);
                w -= hk[k] * v[k];
              }
            hk[i+1] = norm(w);
            if(fabs(hk[i+1]) < 1e-7) // a good breakdown
              v[i+1] = w;
            else
              v[i+1] = w / hk[i+1];  //TODO: how to address hk[i+1] == 0 case?

            {//Given rotations
              for(size_t k = 0; k < i; ++k)
                apply_plane_rotation(hk[k], hk[k+1], cs[k], sn[k]);

              generate_plane_rotation(hk[i],hk[i+1],cs[i],sn[i]);
              apply_plane_rotation(hk[i],hk[i+1],cs[i],sn[i]);
              apply_plane_rotation(s[i],s[i+1],cs[i],sn[i]);

              H[i]=hk;
            }

            if((resid = fabs(s[i+1]/normb)) < eps_){
                update(x_m,i,H,s,v);
                cerr << "# iter " << j << " " << resid << endl;
                cerr << "# converge." << endl;
                return 0;
              }
            cerr << "# iter " << j << " " << resid << endl;
          }

        update(x_m,i-1,H,s,v);
        //r = b_ - A_ * x_m;
        pre->get_r(x, &r[0]);
        rho = norm(r);
        if(fabs(rho/normb) < eps_){
            cerr << "# restart " << rho << endl;
            cerr << "# converge." << endl;
            return 0;
          }
      }
    cerr << "# exceed max iterations." << endl;
    return 1;
  }

  int gmres_solver_sparse::solve(double *x, size_t m, size_t max_iter) const
  {
    itr_matrix<double*> x_m(b_.size(),1,x);
    double resid = 0;
    size_t j = 0;
    matrixd s = zeros<double>(m+1,1);
    matrixd cs(m+1,1), sn(m+1,1);
    matrixd w,r = -1*b_ ;
    hj::sparse::mv(false, A_, x_m, r);
    r *= -1;
    double rho = norm(r);
    s[0] = rho;

    vector<matrixd> v(m+1);
    matrixd hk(m+1,1);
    vector<matrixd> H;
    H.reserve(m+1);
    size_t i;
    while(j < max_iter){
        v[0] = r/rho;
        s *= 0;
        s[0] = rho;

        for(i = 0; i < m && j <= max_iter; ++i, ++j){
            w = zeros<double>(A_.size(1),1);
            hj::sparse::mv(false, A_, v[i], w);
            hk *=0; // clear hk
            for(size_t k = 0; k <= i;++k){
                hk[k] = dot(w, v[k]);
                w -= hk[k] * v[k];
              }
            hk[i+1] = norm(w);
            if(fabs(hk[i+1])<1e-7)
              v[i+1] = w;
            else
              v[i+1] = w / hk[i+1];  //TODO: how to address hk[i+1] == 0 case?

            {//Given rotations
              for(size_t k = 0; k < i; ++k)
                apply_plane_rotation(hk[k], hk[k+1], cs[k], sn[k]);

              generate_plane_rotation(hk[i],hk[i+1],cs[i],sn[i]);
              apply_plane_rotation(hk[i],hk[i+1],cs[i],sn[i]);
              apply_plane_rotation(s[i],s[i+1],cs[i],sn[i]);

              H.push_back(hk);
            }

            if((resid = fabs(s[i+1]/norm(b_))) < eps_){
                update(x_m,i,H,s,v);

               // cerr << "# iter " << j << " " << resid << endl;
                cerr << "# converge." << endl;
                return 0;
              }
            //cerr << "# iter " << j << " " << resid << endl;
          }

        update(x_m,i-1,H,s,v);

        r = -1*b_ ;
        hj::sparse::mv(false, A_, x_m, r);
        r *= -1;

        rho = norm(r);
        if(fabs(rho/norm(b_)) < eps_){
            cerr << "# restart " << rho << endl;
            cerr << "# converge." << endl;
            return 0;
          }
        H.clear();
      }

    cerr << "# exceed max iterations." << endl;
    return 1;
  }

  int gmres_solver_sparse::solve2(
      double *x, size_t m, size_t max_iter,
      const std::shared_ptr<jacobi_precondition_sparse> pre)const
  {
    itr_matrix<double*> x_m(b_.size(),1,x);
    double resid = 0;
    size_t j = 0;
    matrixd s = zeros<double>(m+1,1);
    matrixd cs(m+1,1), sn(m+1,1);
    matrixd w = b_,r = b_;

    matrixd mb = b_;
    pre->get_b(&mb[0]);
    double normb = norm(mb);
    pre->get_r(x, &r[0]); // M*(b-Ax)

    double rho = norm(r);
    s[0] = rho;

    if(fabs(normb) < 1e-7)
      normb = 1.0;

    if((resid = rho/normb) < eps_){
        cerr << "# converge." << endl;
        return 0;
      }

    vector<matrixd> v(m+1);
    matrixd hk(m+1,1);
    vector<matrixd> H(m+1);
    size_t i;
    while(j <= max_iter){
        v[0] = r/rho;
        s *= 0;
        s[0] = rho;

        for(i = 0; i < m && j <= max_iter; ++i, ++j){
            //w = A_ * v[i];
            pre->get_Ax(&v[i][0], &w[0]);
            hk *=0; // clear hk
            for(size_t k = 0; k <= i;++k){
                hk[k] = dot(w, v[k]);
                w -= hk[k] * v[k];
              }
            hk[i+1] = norm(w);
            if(fabs(hk[i+1]) < 1e-7) // a good breakdown
              v[i+1] = w;
            else
              v[i+1] = w / hk[i+1];  //TODO: how to address hk[i+1] == 0 case?

            {//Given rotations
              for(size_t k = 0; k < i; ++k)
                apply_plane_rotation(hk[k], hk[k+1], cs[k], sn[k]);

              generate_plane_rotation(hk[i],hk[i+1],cs[i],sn[i]);
              apply_plane_rotation(hk[i],hk[i+1],cs[i],sn[i]);
              apply_plane_rotation(s[i],s[i+1],cs[i],sn[i]);

              H[i]=hk;
            }

            if((resid = fabs(s[i+1]/normb)) < eps_){
                update(x_m,i,H,s,v);
                cerr << "# iter " << j << " " << resid << endl;
                cerr << "# converge." << endl;
                return 0;
              }
            cerr << "# iter " << j << " " << resid << endl;
          }

        update(x_m,i-1,H,s,v);
        //r = b_ - A_ * x_m;
        pre->get_r(x, &r[0]);
        rho = norm(r);
        if(fabs(rho/normb) < eps_){
            cerr << "# restart " << rho << endl;
            cerr << "# converge." << endl;
            return 0;
          }
      }
    cerr << "# exceed max iterations." << endl;
    return 1;
  }

}
