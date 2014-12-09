#include <iostream>

#include "SQP.h"
#include "qmr.h"

#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/ptree/ptree.h>
#include <zjucad/matrix/io.h>
#include <jtflib/function/function.h>

using namespace std;
using namespace zjucad::matrix;
using namespace hj::sparse;

namespace jtf {

  SQP::SQP()
    :f_(0), c_(0)
  {
    slv_.reset(0);
  }

  int SQP::set_f(jtf::function::functionN1_t<double,int32_t> &f)
  {
    f_ = &f;
    g_ = zeros<double>(f_->dim(), 1);
    s_ = zeros<double>(f_->dim(), 1);
    D_ = zeros<double>(f_->dim(), 1);
    tmp_ = zeros<double>(f_->dim(), 1);
    // do not clean up H if H does not change
    //  H_.resize(0, 0, 0);
    //  corrected_H_.resize(0, 0, 0);
  }

  int SQP::set_c(jtf::function::functionN1_t<double,int32_t> *c)
  {
    c_ = c;
  }

  static double adjust_mu(double mu, double ratio)
  {
    const double min_max_mu[] = {1e-5, 1e6};
    if(ratio < 0) {
        mu *= 4;
        if(mu < min_max_mu[0])
          mu = min_max_mu[0];
      }
    else if(ratio < 0.25) {
        mu *= 4;
        if(mu < min_max_mu[0])
          mu = min_max_mu[0];
      }
    else
      mu /= sqrt(ratio*4); // better than /2
    if(mu > min_max_mu[1])
      mu = min_max_mu[1];
    return mu;
  }

  template <typename VAL_TYPE, typename INT_TYPE>
  inline void add_to_diag(hj::sparse::csc<VAL_TYPE, INT_TYPE> &A, VAL_TYPE *diag)
  {
    for(INT_TYPE ci = 0; ci < A.size(2); ++ci) {
        INT_TYPE nzi = A.ptr()[ci];
        for(; nzi < A.ptr()[ci+1]; ++nzi) {
            if(A.idx()[nzi] == ci) {
                A.val()[nzi] += diag[ci];
                break;
              }
          }
        if(nzi == A.ptr()[ci+1]) {
            cout << "bad add to diag." << endl;
          }
      }
  }

  static void compute_D(double *D, const csc<double, int32_t> &H)
  {
    matrix<double> D1 = zeros<double>(H.size(1), 1);
    for(size_t ci = 0; ci < H.size(2); ++ci) {
        for(size_t nzi = H.ptr()[ci]; nzi < H.ptr()[ci+1]; ++nzi) {
            if(H.idx()[nzi] == ci) {
                D1[ci] = H.val()[nzi];
                if(D1[ci] < 0) {

                    //cerr << "negative H diag: " << D1[ci] << endl;
                    D1[ci] = 0;
                  }
              }
          }
      }
    D1 = sqrt(D1); // STRANGE: should be squared, DTD, according to the
    // book, but removing the sqrt leads to slow
    // convergence.
    for(ptrdiff_t i = 0; i < H.size(1); ++i) {
        if(D[i] < D1[i])
          D[i] = D1[i];
      }
  }

  static void compute_D_jtf(double *D, const csc<double, int32_t> &H)
  {
//    itr_matrix<double*> D0(H.size(1),1, D);
//    D0 = zjucad::matrix::ones<double>(H.size(1), 1);

    matrix<double> D1 = zeros<double>(H.size(1), 1);
    for(size_t ci = 0; ci < H.size(2); ++ci) {
        for(size_t nzi = H.ptr()[ci]; nzi < H.ptr()[ci+1]; ++nzi) {
            if(H.idx()[nzi] == ci) {
                D1[ci] = H.val()[nzi];
                if(D1[ci] < 0) {
                    //cerr << "negative H diag: " << D1[ci] << endl;
                    D1[ci] = -1*D1[ci];
                  }else if(std::fabs(D1[ci]) < 1e-6){
                    D1[ci] = 1;
                  }
              }
          }
      }
    D1 = sqrt(D1); // STRANGE: should be squared, DTD, according to the
    // book, but removing the sqrt leads to slow
    // convergence.
    for(ptrdiff_t i = 0; i < H.size(1); ++i) {
        if(D[i] < D1[i])
          D[i] = D1[i];
      }
  }

  int SQP::solve(double *x, boost::property_tree::ptree &pt)
  {
    if(c_) {
        return solve_by_xx_decomposition(x, pt);
      }

    const double eps[3] = {
      pt.get<double>("epsg.value", 1e-6),
      pt.get<double>("epsf.value", 1e-6),
      pt.get<double>("epsx.value", 1e-20)
    };
    double norm2_of_f_g[2], mu = 1e-6, radius= 1e4;
    pt.put("GS_LM_mu.desc","mu for trust region.");
    if(zjucad::has("GS_LM_mu.value", pt))
      mu = pt.get<double>("GS_LM_mu.value");

    const size_t dim = f_->dim();
    itr_matrix<double *> x0(dim, 1, x);
    int iter_num = pt.get<int>("iter.value");

    for(int i = 0; i < iter_num; ++i) {
        if(cb_)
          if(cb_->at_point(x))
            return __LINE__;

        cerr << endl << i;
        norm2_of_f_g[0] = 0;
        f_->val(x, norm2_of_f_g[0]);
        g_(colon()) = 0;
        f_->gra(x, &g_[0]);
        norm2_of_f_g[1] = dot(g_, g_);

        cerr << "\t" << norm2_of_f_g[0] << "\t" << norm2_of_f_g[1] << "\t" << mu << "\t" << radius;
        if(norm2_of_f_g[1] < eps[0]) {
            cerr << "\ngradient converge." << endl;
            break;
          }

        size_t hes_nnz = 0;
        size_t format = -1;

        if(nnz(H_) == 0 ) { // the first run
            //if(nnz(H_) == 0 || H_.size(1) != dim) { // the first run, JTF modify.
            f_->hes(x, hes_nnz, format, 0, 0, 0);
            H_.resize(dim, dim, hes_nnz);
            //      cerr << "dim: " << f_->dim() << ", nnz: " << hes_nnz << endl;
            f_->hes(x, hes_nnz, format, 0, &H_.ptr()[0], &H_.idx()[0]);
          }
        else
          hes_nnz = nnz(H_);
        H_.val()(colon()) = 0;
        f_->hes(x, hes_nnz, format, &H_.val()[0], &H_.ptr()[0], &H_.idx()[0]);

        bool at_border = true;
        int step_type = -1; // 0 TC, 1 NP, 2 DL
        // Cauchy point
        tmp_ = zeros<double>(dim, 1);
        mv(false, H_, g_, tmp_); // corrected_H?
        const double Cauchy_len = norm2_of_f_g[1]/dot(tmp_, g_), norm_g = sqrt(norm2_of_f_g[1]);
        if(Cauchy_len*norm_g > radius) {
            s_ = -g_*radius/norm_g; // truncate
            step_type = 0;
            cerr << "\tTC";
          }
        else { // Cauchy point is inside, evaluate Newton point
            compute_D(&D_[0], H_);
            D_ *= mu;
            //std::cerr << D_ << std::endl;
            corrected_H_ = H_;
            add_to_diag(corrected_H_, &D_[0]);

            //cerr << "create solver beg" << endl;
            linear_solver * a = linear_solver::create(&corrected_H_.val()[0],
																											&corrected_H_.idx()[0],
																											&corrected_H_.ptr()[0], 
																											corrected_H_.nnz(),
																											corrected_H_.size(1),
																											corrected_H_.size(2),pt);
            slv_.reset(a);
            //cerr << "create solver end" << endl;
            if(slv_.get()) {
                slv_->solve(&g_[0], &s_[0], 1, pt);
                s_ = -s_;

                if(norm(s_) < radius) { // Newton point is inside of the radius
                    //        x0 += s_;
                    step_type = 1;
                    cerr << "\tNP";
                    at_border = false;
                  }
                else { // intersect with radius
                    const matrix<double> CP = -g_*Cauchy_len;
                    const double gamma = dot(CP, CP)/dot(s_, s_);
                    matrix<double> C2N = (0.8*gamma+0.2)*s_-CP;
                    const double quadric_eq[3] = {dot(C2N, C2N), 2*dot(CP, C2N), dot(CP, CP)-radius*radius};
                    double Delta = quadric_eq[1]*quadric_eq[1]-4*quadric_eq[0]*quadric_eq[2];
                    assert(Delta > 0);
                    Delta = sqrt(Delta);
                    if(quadric_eq[0] < 0) {
                        cerr << "\tbad dog leg: " << quadric_eq[0] << endl;
                      }
                    const double roots[2] = {
                      (-quadric_eq[1]+Delta)/(2*quadric_eq[0]),
                      (-quadric_eq[1]-Delta)/(2*quadric_eq[0])
                    };
                    // it's highly possible that dot(CP, C2N) > 0
                    if(roots[0] >= 0 && roots[0] <= 1) {
                        if(roots[1] >= 0 && roots[1] <= 1)
                          cerr << "two intersection?" << endl;
                        s_ = CP+roots[0]*C2N;
                      }
                    else if(roots[1] >= 0 && roots[1] <= 1) {
                        cerr << "strange why the second root?" << endl;
                        s_ = CP+roots[1]*C2N;
                      }
                    else {
                        if(roots[0] > 0) {
                            s_ = CP+roots[0]*C2N;
                          }
                        if(roots[1] > 0) {
                            s_ = CP+roots[1]*C2N;
                          }
                        if(roots[0]*roots[1] > 0)
                          cerr << "\nunbelievable that no root is OK: " << roots[0] << " " << roots[1] << endl;
                      }
                    step_type = 2;
                    cerr << "\tDL";
                  }
              }
            else { // H is not SPD
                mu *= 16;
                continue;
              }
          }
        matrix<double> ori_x = x0;
        if(max(fabs(s_)) < eps[2]) {
            cerr << "\nstep converge." << endl;
            break;
          }
        x0 += s_;
        // update mu and radius
        double new_val = 0, est_val;
        f_->val(x, new_val);
        tmp_(colon()) = 0;
        mv(false, H_, s_, tmp_);
        est_val = norm2_of_f_g[0] + dot(s_, g_) + dot(s_, tmp_)/2;
        double ratio = -1;
        if(new_val < norm2_of_f_g[0]) {
            const double delta_f = new_val-norm2_of_f_g[0];
            //			if(delta_f < 0 && -delta_f < eps[1]) {
            //        cerr << "\nx converge." << endl;
            //				break;
            //      }
            double t = est_val-norm2_of_f_g[0];
            if(fabs(t) < 1e-10)
              t = 1e-5;
            ratio = delta_f/t;
          }
        //    cerr << "\t" << est_val << "\t" << new_val << "\t" << norm2_of_f_g[0] << "\t" << ratio;
        if(ratio < 0 && step_type == 0) { // accept bad step for NP and DL
            x0 = ori_x;
          }
        mu = adjust_mu(mu, ratio);
        if(ratio < 0.25) {
            radius /= 4;
          }
        else if(ratio > 0.75 || at_border) {
            radius *= 2;
          }
        cerr << endl;
      }
    cerr << "\noptimization over." << endl;
    return 0;
  }

  int SQP::solve_by_xx_decomposition(double *x, boost::property_tree::ptree &pt)
  {
    const double eps[3] = {
      pt.get<double>("epsg.value", 1e-6),
      pt.get<double>("epsf.value", 1e-6),
      pt.get<double>("epsx.value", 1e-20)
    };
    double norm2_of_f_g[2], mu = 1e-6, radius= 1e4;
    pt.put("GS_LM_mu.desc","mu for trust region.");
    if(zjucad::has("GS_LM_mu.value", pt))
      mu = pt.get<double>("GS_LM_mu.value");

    const size_t dim = f_->dim();
    itr_matrix<double *> x0(dim, 1, x);
    tmp_ = zeros<double>(dim, 1);
    matrix<double> nabla_c(dim, 1);

    int iter_num = pt.get<int>("iter.value");

    double lambda = 0;
    for(int i = 0; i < iter_num; ++i) {
        if(cb_)
          if(cb_->at_point(x))
            return __LINE__;
        cerr << i;

        norm2_of_f_g[0] = 0;
        f_->val(x, norm2_of_f_g[0]);
        g_(colon()) = 0;
        f_->gra(x, &g_[0]);
        norm2_of_f_g[1] = dot(g_, g_);

        double c = 0;
        c_->val(x, c);
        nabla_c(colon()) = 0;
        c_->gra(x, &nabla_c[0]);

        const double total_gradient = std::max(max(fabs(g_ + lambda * nabla_c)), fabs(c));
        cerr << "\tf,g,g+c: " << norm2_of_f_g[0] << "\t" << norm2_of_f_g[1] << "\t" << total_gradient << "\t" << mu << "\t" << radius;
        if(total_gradient < eps[0]) {
            cerr << "\ngradient converge." << total_gradient << " " << eps[0] << endl;
            break;
          }

        size_t hes_nnz = 0, format = -1;
        if(nnz(H_) == 0) { // the first run
            clock_t beg = clock();
            f_->hes(x, hes_nnz, format, 0, 0, 0);
            cout << "\n# hes query: " << (clock()-beg)/double(CLOCKS_PER_SEC) << endl;
            H_.resize(dim, dim, hes_nnz);
            //      cerr << "dim: " << f_->dim() << ", nnz: " << hes_nnz << endl;
            f_->hes(x, hes_nnz, format, 0, &H_.ptr()[0], &H_.idx()[0]);
          }
        else
          hes_nnz = nnz(H_);
        H_.val()(colon()) = 0;
        f_->hes(x, hes_nnz, format, &H_.val()[0], &H_.ptr()[0], &H_.idx()[0]);
        compute_D(&D_[0], H_);

        D_ *= mu;
        corrected_H_ = H_;
        add_to_diag(corrected_H_, &D_[0]);

        //cerr << "create solver beg" << endl;
        slv_.reset(linear_solver::create(&corrected_H_.val()[0], &corrected_H_.idx()[0],
																				 &corrected_H_.ptr()[0], corrected_H_.nnz(), 
																				 corrected_H_.size(1), corrected_H_.size(2),pt));
        //cerr << "create solver end" << endl;
        if(slv_.get()) {
            //      cerr << "\nconstraint: " << c << " " << norm(nabla_c) << endl;
            slv_->solve(&g_[0], &s_[0], 1, pt);
            slv_->solve(&nabla_c[0], &tmp_[0], 1, pt);
            double xx = dot(tmp_, nabla_c);
            if(xx < 0) {
                cerr << "# [info] xx " << xx << endl;
                cerr << "impossible" << endl;
                //exit(0);
              }
            // if(xx < 1e-10) {
            //   cerr << "\ninstability caused by g_c*inv(H)*g_c: " << c << " " << norm(nabla_c) << endl;
            //   break;
            //   xx = 1e10;
            // }
            lambda = (c - dot(s_, nabla_c))/xx;
            if(!isfinite(lambda)) {
                cerr << "nan or inf stop" << endl;
                return 0;
              }
            cout << "\tlambda: " << lambda << endl;
            s_ = -s_ - lambda*tmp_;
          }
        else { // H is not SPD
            mu *= 16;
            continue;
          }

        matrix<double> ori_x = x0;
        if(max(fabs(s_)) < eps[2]) {
            cerr << "\nstep converge." << endl;
            break;
          }
        const double s_len = norm(s_);
        if(s_len > radius) {
            s_ *= radius/s_len;
            cout << "clamp: " << s_len << endl;
          }
        size_t k, MAX_INTERIOR = 4;
        for(k = 0; k < MAX_INTERIOR; ++k) {
            x0 = ori_x + s_;
            if(f_->is_valid(x)) break;
            s_ /= 2;
            cerr << "# scale for interior: " << k << "\r" << flush;
          }
        if(k == MAX_INTERIOR) {
            cerr << "# interior fail." << endl;
            return -1; // jtf: if interiro fail, neet to remesh topology or exit
          }
        // update mu and radius
        double new_val = 0, est_val;
        f_->val(x, new_val);
        tmp_(colon()) = 0;
        mv(false, H_, s_, tmp_);
        est_val = norm2_of_f_g[0] + dot(s_, g_) + dot(s_, tmp_)/2;
        double ratio = -1;
        if(new_val < norm2_of_f_g[0]) {
            const double delta_f = new_val-norm2_of_f_g[0];
            double t = est_val-norm2_of_f_g[0];
            if(fabs(t) < 1e-10)
              t = 1e-5;
            ratio = delta_f/t;
          }
        //    cerr << "\t" << est_val << "\t" << new_val << "\t" << norm2_of_f_g[0] << "\t" << ratio;
        if(ratio < 0) {
            x0 = ori_x;
            radius = norm(s_);
            const double dir = dot(g_, s_)/(norm(g_)*radius);
            cout << "bad step, direction: " << dir << endl;
            if(dir > 0 && mu < 10)
              mu = 10;
          }
        mu = adjust_mu(mu, ratio);
        if(ratio < 0.25) {
            radius /= 4;
          }
        else if(ratio > 0.75) {
            radius *= 2;
          }
        cerr << endl;
      }
    cerr << "\noptimization over." << endl;
    return 0;
  }

  SQPEX::SQPEX()
  {

  }

  /// \brief SQPEX::solve
  /// \param x
  /// \param pt
  /// \return
  int SQPEX::solve(double *x, boost::property_tree::ptree &pt)
  {
    H_.resize(0,0);

    if(c_) {
        return solve_by_xx_decomposition(x, pt);
      }

    const double eps[3] = {
      pt.get<double>("epsg.value", 1e-8),
      pt.get<double>("epsf.value", 1e-8),
      pt.get<double>("epsx.value", 1e-8)
    };

    const size_t dim = f_->dim();
    itr_matrix<double *> x0(dim, 1, x);
    double mu = 1e-3;
    pt.put("GS_LM_mu.desc","mu for trust region.");
    if(zjucad::has("GS_LM_mu.value", pt)){
        mu = pt.get<double>("GS_LM_mu.value");
      }

    double norm2_of_f_g[2];

    pt.put("dump_step.desc","dump every n steps");
    size_t dump_step = pt.get<size_t>("dump_step.value", size_t(0));
    string dump_pref;
    pt.put("dump_pref.desc","prefix of the files for dumping x");
    if(dump_step > 0)
      dump_pref = pt.get<string>("dump_pref.value");

    double sigma = 1.0,cur_mu = 0;
    const int iter_num = pt.get<int>("iter.value");
    const double min_mu = 1e-10,max_mu = 1e8,min_sigma = 1e-20;

    int i,result = 0;
    size_t format;
    for(i = 0; i < iter_num; ++i) {
        cout<<"inner iteration : "<<i+1<<std::endl;

        norm2_of_f_g[0] = 0;
        f_->val(x, norm2_of_f_g[0]);
        g_(colon()) = 0;
        f_->gra(x, &g_[0]);
        norm2_of_f_g[1] = dot(g_, g_);

        cerr << "\t" << norm2_of_f_g[0] << "\t" << norm2_of_f_g[1] << "\t" << mu << "\n";
        if(norm2_of_f_g[1] < eps[0]){
            std::cout<<"gradient convergence!"<<std::endl;
            result = 1;
            break;
          }

        //cout << "mu: " << mu << endl;
        size_t hes_nnz = 0;
        if(nnz(H_) == 0 ) { // the first run
            //if(nnz(H_) == 0 || H_.size(1) != dim) { // the first run, JTF modify.

            f_->hes(x, hes_nnz, format, 0, 0, 0);
            H_.resize(dim, dim, hes_nnz);
            //      cerr << "dim: " << f_->dim() << ", nnz: " << hes_nnz << endl;
            f_->hes(x, hes_nnz,format, 0, &H_.ptr()[0], &H_.idx()[0]);
          }
        else
          hes_nnz = nnz(H_);
        H_.val()(colon()) = 0;
        f_->hes(x, hes_nnz, format, &H_.val()[0], &H_.ptr()[0], &H_.idx()[0]);

        bool success = true;
        matrix<double> ori_x = x0;
        while(mu < max_mu){
            compute_D(&D_[0], H_);
            D_ *= mu;
            corrected_H_ = H_;
            add_to_diag(corrected_H_, &D_[0]);

            slv_.reset(linear_solver::create(&corrected_H_.val()[0],
																						 &corrected_H_.idx()[0],
																						 &corrected_H_.ptr()[0],
																						 corrected_H_.nnz(),
																						 corrected_H_.size(1),
																						 corrected_H_.size(2), pt));

            if(slv_.get()) {
                slv_->solve(&g_[0], &s_[0], 1, pt);
                s_ = -s_;

                if(mu > min_mu){
                    mu *= 0.1;
                  }

                double new_val(0);
                x0 = ori_x + s_;
                f_->val(x, new_val);
                tmp_(colon()) = 0;
                mv(false, H_, s_, tmp_);
                double est_val = norm2_of_f_g[0] + dot(s_, g_) + dot(s_, tmp_)/2;
                double ratio = -1;
                if(new_val < norm2_of_f_g[0]) {
                    const double delta_f = new_val-norm2_of_f_g[0];
                    double t = est_val-norm2_of_f_g[0];
                    if(fabs(t) < 1e-10)
                      t = 1e-5;
                    ratio = delta_f/t;
                  }

                if(ratio < 0.25) {
                    mu *= 3;
                    if(mu < min_mu)
                      mu = min_mu;
                  }
                else{
                    mu /= 2*sqrt(ratio); // better than /2
                  }
                if(mu > max_mu)
                  mu = max_mu;
                cur_mu = success ? 0 : mu;

                break;
              }else{
                if(mu < max_mu){
                    mu *= 10;
                  }
                success = false;
              }
          }

        if(max(fabs(s_)) < eps[2]) {
            cerr << "\nstep converge." << endl;
            result = 0;
            break;
          }

        double new_val;
        sigma = std::min(sigma,1.0);
        size_t tried_num = 0;
        const size_t max_tried_num = 4;
        for(; tried_num < max_tried_num; ++tried_num){
            x0 = ori_x + sigma*s_;
            if(f_->is_valid(&x0[0])) break;
            sigma *= 0.5;
          }

        if(tried_num == max_tried_num){
            cerr << "# interiro fail." << endl;
            return -1;
          }

        while(1) {
            x0 = ori_x + sigma*s_;
            new_val = 0;
            f_->val(x, new_val);
            if(norm2_of_f_g[0] > new_val ||
               sigma < min_sigma){
                sigma *= 2;
                break;
              }
            else{
                sigma *= 0.5;
                tried_num++;
              }
          }
      }

    if(iter_num == i){
        result = -2;
      }

    std::cerr << "\noptimization over." << endl;

    return result;
  }


  int SQPEX::solve_by_xx_decomposition(double *x, boost::property_tree::ptree &pt)
  {
    const double eps[3] = {
      pt.get<double>("epsg.value", 1e-6),
      pt.get<double>("epsf.value", 1e-6),
      pt.get<double>("epsx.value", 1e-20)
    };
    double norm2_of_f_g[2], mu = 1e-6, radius= 1e4;
    pt.put("GS_LM_mu.desc","mu for trust region.");
    if(zjucad::has("GS_LM_mu.value", pt))
      mu = pt.get<double>("GS_LM_mu.value");

    const size_t dim = f_->dim();
    itr_matrix<double *> x0(dim, 1, x);
    tmp_ = zeros<double>(dim, 1);
    matrix<double> nabla_c(dim, 1);

    double sigma = 1.0,cur_mu = 0,lambda = 0;
    const int iter_num = pt.get<int>("iter.value");
    const double min_mu = 1e-10,max_mu = 1e8,min_sigma = 1e-20;
    size_t format;
    for(int i = 0; i < iter_num; ++i) {
        if(cb_)
          if(cb_->at_point(x))
            return __LINE__;
        cerr << i;

        norm2_of_f_g[0] = 0;
        f_->val(x, norm2_of_f_g[0]);
        g_(colon()) = 0;
        f_->gra(x, &g_[0]);
        norm2_of_f_g[1] = dot(g_, g_);

        double c = 0;
        c_->val(x, c);
        nabla_c(colon()) = 0;
        c_->gra(x, &nabla_c[0]);

        const double total_gradient = std::max(max(fabs(g_ + lambda * nabla_c)), fabs(c));
        cerr << "\tf,g,g+c: " << norm2_of_f_g[0] << "\t" <<
                norm2_of_f_g[1] << "\t" << total_gradient << "\t" << mu << "\t" << radius;
        if(total_gradient < eps[0]) {
            cerr << "\ngradient converge." << total_gradient << " " << eps[0] << endl;
            break;
          }

        size_t hes_nnz = 0;
        if(nnz(H_) == 0) { // the first run
            f_->hes(x, hes_nnz, format, 0, 0, 0);
            H_.resize(dim, dim, hes_nnz);
            //      cerr << "dim: " << f_->dim() << ", nnz: " << hes_nnz << endl;
            f_->hes(x, hes_nnz, format, 0, &H_.ptr()[0], &H_.idx()[0]);
          }
        else
          hes_nnz = nnz(H_);
        H_.val()(colon()) = 0;
        f_->hes(x, hes_nnz,format, &H_.val()[0], &H_.ptr()[0], &H_.idx()[0]);

        bool success = true;
        matrix<double> ori_x = x0;
        while(mu < max_mu){
            compute_D(&D_[0], H_);
            D_ *= mu;
            corrected_H_ = H_;
            add_to_diag(corrected_H_, &D_[0]);

            slv_.reset(linear_solver::create(&corrected_H_.val()[0], 
																						 &corrected_H_.idx()[0], 
																						 &corrected_H_.ptr()[0],
																						 corrected_H_.nnz(),
																						 corrected_H_.size(1),
																						 corrected_H_.size(2),pt));
            //cerr << "create solver end" << endl;
            if(slv_.get()) {
                //      cerr << "\nconstraint: " << c << " " << norm(nabla_c) << endl;
                slv_->solve(&g_[0], &s_[0], 1, pt);
                slv_->solve(&nabla_c[0], &tmp_[0], 1, pt);
                double xx = dot(tmp_, nabla_c);
                if(xx < 0) {
                    cerr << "# [info] xx " << xx << endl;
                    cerr << "impossible" << endl;
                  }

                lambda = (c - dot(s_, nabla_c))/xx;
                if(!isfinite(lambda)) {
                    cerr << "nan or inf stop" << endl;
                    return 0;
                  }
                cout << "\tlambda: " << lambda << endl;
                s_ = -s_ - lambda*tmp_;

                if(mu > min_mu){
                    mu *= 0.1;
                  }

                double new_val(0);
                x0 = ori_x + s_;
                f_->val(x, new_val);
                tmp_(colon()) = 0;
                mv(false, H_, s_, tmp_);
                double est_val = norm2_of_f_g[0] + dot(s_, g_) + dot(s_, tmp_)/2;
                double ratio = -1;
                if(new_val < norm2_of_f_g[0]) {
                    const double delta_f = new_val-norm2_of_f_g[0];
                    double t = est_val-norm2_of_f_g[0];
                    if(fabs(t) < 1e-10)
                      t = 1e-5;
                    ratio = delta_f/t;
                  }

                if(ratio < 0.25) {
                    mu *= 3;
                    if(mu < min_mu)
                      mu = min_mu;
                  }
                else{
                    mu /= 2*sqrt(ratio); // better than /2
                  }
                if(mu > max_mu)
                  mu = max_mu;

                cur_mu = success ? 0 : mu;

                break;
              }else{
                if(mu < max_mu){
                    mu *= 10;
                  }
                success = false;
              }
          }

        if(max(fabs(s_)) < eps[2]) {
            cerr << "\nstep converge." << endl;
            break;
          }

        double new_val;
        sigma = std::min(sigma,1.0);
        size_t tried_num = 0;
        const size_t max_tried_num = 4;
        for(; tried_num < max_tried_num; ++tried_num){
            x0 = ori_x + sigma*s_;
            if(f_->is_valid(&x0[0])) break;
            sigma *= 0.5;
          }

        if(tried_num == max_tried_num){
            cerr << "# interiro fail." << endl;
            return -1;
          }

        while(1) {
            x0 = ori_x + sigma*s_;
            new_val = 0;
            f_->val(x, new_val);
            if(norm2_of_f_g[0] > new_val ||
               sigma < min_sigma){
                sigma *= 2;
                break;
              }else{
                sigma *= 0.5;
              }
          }

        cerr << endl;
      }
    cerr << "\noptimization over." << endl;
    return 0;
  }

}
