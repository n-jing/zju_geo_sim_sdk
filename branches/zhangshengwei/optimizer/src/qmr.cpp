#include "qmr.h"
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <iostream>

using namespace std;
using namespace zjucad::matrix;

namespace jtf{

  int qmr_solver::solve(double *x_raw, size_t max_iter) const
  {
    double resid;
    itr_matrix<double*> x(A_.size(2),1,x_raw);

    double rho, rho_1, xi, gamma, gamma_1, theta, theta_1;
    double eta, delta, ep, beta;
    matrix<double> r, v_tld, y, w_tld, z;
    matrix<double> v, w, y_tld, z_tld;
    matrix<double> p, q, p_tld, d, s;

    double normb = norm(b_);

    r = b_ - A_ * x;

    if (fabs(normb) < 1e-7)
      normb = 1;

    if ((resid = norm(r) / normb) <= eps_) {
        cerr << "# [info] converge." << endl;
        return 0;
      }

    v_tld = r;
    y = v_tld;
    rho = norm(y);

    w_tld = r;
    z = w_tld;
    xi = norm(z);

    gamma = 1.0;
    eta = -1.0;
    theta = 0.0;

    for (int i = 1; i <= max_iter; i++) {

        if (fabs(rho) < 1e-7){
            cerr << "# breakdown: 2" << endl;
            return 2;                        // return on breakdown
          }

        if (fabs(xi) < 1e-7){
            cerr << "# breakdown 7" << endl;
            return 7;                        // return on breakdown
          }

        v = (1. / rho) * v_tld;
        y = (1. / rho) * y;

        w = (1. / xi) * w_tld;
        z = (1. / xi) * z;

        delta = dot(z, y);
        if (fabs(delta) < 1e-7){
            cerr << "# breakdown 5" << endl;
            return 5;                        // return on breakdown
          }

        y_tld = y;               // apply preconditioners
        z_tld = z;

        if (i > 1) {
            p = y_tld - (xi * delta / ep) * p;
            q = z_tld - (rho * delta / ep) * q;
          } else {
            p = y_tld;
            q = z_tld;
          }

        p_tld = A_ * p;
        ep = dot(q, p_tld);
        if (fabs(ep) < 1e-7){
            cerr << "# breakdown 6" << endl;
            return 6;                        // return on breakdown
          }

        beta = ep / delta;
        if (fabs(beta) < 1e-7){
            cerr << "# breakdown 3" << endl;
            return 3;                        // return on breakdown
          }

        v_tld = p_tld - beta * v;
        y = v_tld;

        rho_1 = rho;
        rho = norm(y);
        w_tld = trans(A_) * q - beta * w;
        z = w_tld;

        xi = norm(z);

        gamma_1 = gamma;
        theta_1 = theta;

        theta = rho / (gamma_1 * beta);
        gamma = 1.0 / sqrt(1.0 + theta * theta);

        if (fabs(gamma) < 1e-7){
            cerr << "# breakdown 4" << endl;
            return 4;                        // return on breakdown
          }

        eta = -eta * rho_1 * gamma * gamma /
            (beta * gamma_1 * gamma_1);

        if (i > 1) {
            d = eta * p + (theta_1 * theta_1 * gamma * gamma) * d;
            s = eta * p_tld + (theta_1 * theta_1 * gamma * gamma) * s;
          } else {
            d = eta * p;
            s = eta * p_tld;
          }

        x += d;                            // update approximation vector
        r -= s;                            // compute residual

        if ((resid = norm(r) / normb) <= eps_) {
            cerr << "# [info] converge." << endl;
            return 0;
          }
        cerr << " # iter " << i << " " << resid << endl;
      }
    cerr << "# [info] exceed max iterations."  << endl;
    return 1;                            // no convergence
  }

  int qmr_solver_sparse::solve(double *x_raw, size_t max_iter) const
  {
    double resid;
    itr_matrix<double*> x(A_.size(2),1,x_raw);

    double rho, rho_1, xi, gamma, gamma_1, theta, theta_1;
    double eta, delta, ep, beta;
    matrix<double> r, v_tld, y, w_tld, z;
    matrix<double> v, w, y_tld, z_tld;
    matrix<double> p, q, p_tld, d, s;

    double normb = norm(b_);

    r = -1 *b_;
    hj::sparse::mv(false, A_, x, r);
    r *= -1;

    if (fabs(normb) < 1e-7)
      normb = 1;

    if ((resid = norm(r) / normb) <= eps_) {
        cerr << "# [info] converge. iter 0" << endl;
        return 0;
      }

    v_tld = r;
    y = v_tld;
    rho = norm(y);

    w_tld = r;
    z = w_tld;
    xi = norm(z);

    gamma = 1.0;
    eta = -1.0;
    theta = 0.0;

    for (int i = 1; i <= max_iter; i++) {

        if (fabs(rho) < 1e-7){
            cerr << "# breakdown: 2 iter " << i << endl;
            return 2;                        // return on breakdown
          }

        if (fabs(xi) < 1e-7){
            cerr << "# breakdown 7 iter " << i << endl;
            return 7;                        // return on breakdown
          }

        v = (1. / rho) * v_tld;
        y = (1. / rho) * y;

        w = (1. / xi) * w_tld;
        z = (1. / xi) * z;

        delta = dot(z, y);
        if (fabs(delta) < 1e-7){
            cerr << "# breakdown 5 iter " << i << endl;
            return 5;                        // return on breakdown
          }

        y_tld = y;               // apply preconditioners
        z_tld = z;

        if (i > 1) {
            p = y_tld - (xi * delta / ep) * p;
            q = z_tld - (rho * delta / ep) * q;
          } else {
            p = y_tld;
            q = z_tld;
          }

        p_tld = zeros<double>(p.size(),1);
        hj::sparse::mv(false, A_,p,p_tld);

        ep = dot(q, p_tld);
        if (fabs(ep) < 1e-7){
            cerr << "# breakdown 6 iter " << i << endl;
            return 6;                        // return on breakdown
          }

        beta = ep / delta;
        if (fabs(beta) < 1e-7){
            cerr << "# breakdown 3 iter " << i << endl;
            return 3;                        // return on breakdown
          }

        v_tld = p_tld - beta * v;
        y = v_tld;

        rho_1 = rho;
        rho = norm(y);
        w_tld = -1*beta *w;
        hj::sparse::mv(true, A_, q, w_tld);

        z = w_tld;

        xi = norm(z);

        gamma_1 = gamma;
        theta_1 = theta;

        theta = rho / (gamma_1 * beta);
        gamma = 1.0 / sqrt(1.0 + theta * theta);

        if (fabs(gamma) < 1e-7){
            cerr << "# breakdown 4: iter " << i << endl;
            return 4;                        // return on breakdown
          }

        eta = -eta * rho_1 * gamma * gamma /
            (beta * gamma_1 * gamma_1);

        if (i > 1) {
            d = eta * p + (theta_1 * theta_1 * gamma * gamma) * d;
            s = eta * p_tld + (theta_1 * theta_1 * gamma * gamma) * s;
          } else {
            d = eta * p;
            s = eta * p_tld;
          }

        x += d;                            // update approximation vector
        r -= s;                            // compute residual

        if ((resid = norm(r) / normb) <= eps_) {
            cerr << "# [info] converge. iter " << i << endl;
            return 0;
          }
       // cerr << " # iter " << i << " " << resid << endl;
      }
    cerr << "# [info] exceed max iterations. iter " << max_iter  << endl;
    return 1;                            // no convergence
  }
}
