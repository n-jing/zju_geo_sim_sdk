#ifndef THIRD_SOLVER_IPOPT_H
#define THIRD_SOLVER_IPOPT_H

#include <Ipopt/coin/IpIpoptApplication.hpp>
#include <Ipopt/coin/IpTNLP.hpp>

#include <boost/shared_ptr.hpp>
#include <zjucad/matrix/itr_matrix.h>
#include "kernel_zjucad_matrix.h"
#include "objective_function.h"

namespace zjucad
{
namespace LSCO
{

template <typename MATRIX_TYPE>
struct kernel_traits;

class solver_use_ipopt : public Ipopt::TNLP
{
public:
    typedef kernel_traits<hj::sparse::csc<Ipopt::Number,int32_t> > kernel_type;
    typedef boost::shared_ptr<const objective_function<kernel_type> > func_ptr;

    explicit solver_use_ipopt(func_ptr func,
                              zjucad::matrix::matrix<Ipopt::Number> &init_val)
        :func_(func), init_val_(init_val) {}
    virtual ~solver_use_ipopt() {}

public:
    virtual bool get_nlp_info(
        Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style) {
        n = func_->nr_inputs();
        // contains equality constraint and inequality constraint
        m = func_->nr_equals()+func_->nr_inequals();
        //std::vector<double> x(n);
        static zjucad::matrix::matrix<Ipopt::Number> fake_x(n,1);
        nnz_jac_g = 0;
        kernel_type::matrix_builder mb_jac, mb_hes;
        func_->eval_constraint_gradient(fake_x, mb_jac , false);

        nnz_jac_g = mb_jac.size();

        zjucad::matrix::matrix<Ipopt::Number> lambda(m,1);
        Ipopt::Number coef = 1.0;
        func_->eval_lagrangian_hessian(fake_x, lambda, coef, mb_hes,false);
        nnz_h_lag = mb_hes.size();

        index_style = TNLP::C_STYLE;

        return true;
    }

    //! @brief Method to return the bounds for my problem
    virtual bool get_bounds_info(
        Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
        Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) {
        //check the value of n and m
        assert(n == func_->nr_inputs());
        assert(m == func_->nr_equals() + func_->nr_inequals());

        // The default value of nlp_upper_bound_inf and nlp_lower_bound_inf
        // is 1e19 and can be changed through ipopt options.

        // the variables have none lower or upper bounds
        for (Ipopt::Index i = 0; i < n; ++i) {
            x_l[i] = -1e19;
            x_u[i] = 1e19;
        }

        // equality constraint
        Ipopt::Index i = 0;
        for (; i < func_->nr_equals(); ++i) {
            g_l[i] = 0;
            g_u[i] = 0;
        }

        for(; i < func_->nr_equals() + func_->nr_inequals(); ++i) {
            g_l[i] = 0;
            g_u[i] = 1e19;
        }

        return true;
    }

    //! @brief Method to return the starting point for the algorithm
    virtual bool get_starting_point(
        Ipopt::Index n, bool init_x, Ipopt::Number* x,
        bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
        Ipopt::Index m, bool init_lambda,
        Ipopt::Number* lambda) {
        assert(init_x == true);
        assert(init_z == false);
        assert(init_lambda == false);

        std::copy(init_val_.begin(), init_val_.end(), x);

        return true;
    }

    //! @brief Method to return the objective value
    virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                        Ipopt::Number& obj_value) {
        assert(n == func_->nr_inputs());

        zjucad::matrix::itr_matrix<const Ipopt::Number*> x0(n,1,x);
        obj_value = func_->eval(x0,false);

        return true;
    }

    //! @brief Method to return the gradient of the objective
    virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                             Ipopt::Number* grad_f) {
        assert(n == func_->nr_inputs());
        zjucad::matrix::itr_matrix<Ipopt::Number*> grad_f_m(func_->nr_inputs(),1,grad_f);
        zjucad::matrix::itr_matrix<const Ipopt::Number*> x0(n,1,x);

        zjucad::matrix::matrix<Ipopt::Number> grad_result(func_->nr_inputs(),1);

        func_->eval_func_gradient(x0,grad_result,true);

        grad_f_m = grad_result;

        return true;
    }

    //! @brief Method to return the constraint residuals
    virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                        Ipopt::Index m, Ipopt::Number* g) {
        if(func_->nr_equals() == 0 && func_->nr_inequals() == 0)
            return true;

        assert(n == func_->nr_inputs());
        assert(m == func_->nr_equals() + func_->nr_inequals());

        zjucad::matrix::itr_matrix<const Ipopt::Number*> x0(n,1,x);
        zjucad::matrix::matrix<Ipopt::Number> x_m = x0;
        zjucad::matrix::itr_matrix<Ipopt::Number*> g0(m,1,g);
        zjucad::matrix::matrix<Ipopt::Number> g_m = g0;

        func_->eval_constraint_value(x_m,g_m, new_x);

        g0 = g_m;
        return true;
    }

    //! @brief The structure of the jacobian (if "values" is NULL)
    //!        The values of the jacobian (if "values" is not NULL)
    virtual bool eval_jac_g(
        Ipopt::Index n, const Ipopt::Number* x, bool new_x,
        Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
        Ipopt::Index *jCol, Ipopt::Number* values) {
        kernel_type::matrix_builder mb;
        if(values == NULL) {
            static zjucad::matrix::matrix<Ipopt::Number> fake_x(n,1);
            func_->eval_constraint_gradient(fake_x,mb,false);
            for(size_t i = 0; i < mb.size(); ++i) {
                iRow[i] = mb[i].first.first;
                jCol[i] = mb[i].first.second;
            }
        } else {
            zjucad::matrix::itr_matrix<const Ipopt::Number*> x0(n,1,x);
            zjucad::matrix::matrix<Ipopt::Number> fake_x = x0;
            func_->eval_constraint_gradient(fake_x,mb,false);
            for(size_t i = 0; i < mb.size(); ++i) {
                values[i] = mb[i].second;
            }
        }
        return true;
    }

    //! @brief The structure of the hessian of the lagrangian (if "values" is NULL)
    //!        The values of the hessian of the lagrangian (if "values" is not NULL)
    virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                        Ipopt::Number obj_factor, Ipopt::Index m,
                        const Ipopt::Number* lambda, bool new_lambda,
                        Ipopt::Index nele_hess, Ipopt::Index* iRow,
                        Ipopt::Index* jCol, Ipopt::Number* values) {
        kernel_type::matrix_builder hessian;
        if(values == NULL) {
            static zjucad::matrix::matrix<Ipopt::Number> fake_x(n,1);
            zjucad::matrix::matrix<Ipopt::Number> lambda_m(m,1);
            func_->eval_lagrangian_hessian(fake_x, lambda_m, obj_factor, hessian, false);
            for(size_t i = 0; i < hessian.size(); ++i) {
                iRow[i] = hessian[i].first.first;
                jCol[i] = hessian[i].first.second;
            }
        } else {
            zjucad::matrix::itr_matrix<const Ipopt::Number*> x0(n,1,x);
            zjucad::matrix::matrix<Ipopt::Number> fake_x = x0;
            zjucad::matrix::itr_matrix<const Ipopt::Number*> lambda0(m,1,lambda);
            zjucad::matrix::matrix<Ipopt::Number> fake_lambda = lambda0;
            func_->eval_lagrangian_hessian(fake_x, fake_lambda, obj_factor, hessian, false);
            for(size_t i = 0; i < hessian.size(); ++i) {
                values[i] = hessian[i].second;
            }
        }

        return true;
    }

    //! @brief  This method is called when the algorithm is complete
    //!         so the TNLP can store/write the solution */
    virtual void finalize_solution(
        Ipopt::SolverReturn status,
        Ipopt::Index n, const Ipopt::Number* x,
        const Ipopt::Number* z_L, const Ipopt::Number* z_U,
        Ipopt::Index m, const Ipopt::Number* g,
        const Ipopt::Number* lambda,
        Ipopt::Number obj_value,
        const Ipopt::IpoptData* ip_data,
        Ipopt::IpoptCalculatedQuantities* ip_cq) {
        assert(n == init_val_.size());
        std::copy(x, x+ n, init_val_.begin());
    }

private:
    func_ptr func_;
    zjucad::matrix::matrix<Ipopt::Number> &init_val_;
private:
    solver_use_ipopt(const solver_use_ipopt&);

    /** Overloaded Equals Operator */
    void operator=(const solver_use_ipopt&);
};

class third_solver_ipopt
{
public:
    typedef kernel_traits<hj::sparse::csc<Ipopt::Number,int32_t> > kernel_type;
    typedef boost::shared_ptr<objective_function<kernel_type> > func_ptr;

    third_solver_ipopt(func_ptr & obj):obj_(obj) {}

    int solve(zjucad::matrix::matrix<double> & init_x,
              boost::property_tree::ptree & pt) {
        const size_t max_iter_num = pt.get<size_t>("iter.value");
        Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
        app->Options()->SetNumericValue("tol", 1e-7);
        //app->Options()->SetStringValue("derivative_test", "second-order");
        app->Options()->SetStringValue("output_file", "ipopt.out");
        //app->Options()->SetStringValue("mu_strategy", "adaptive");
        app->Options()->SetIntegerValue("max_iter", max_iter_num);
        app->Options()->SetStringValue("linear_solver", "ma97");
        app->Options()->SetNumericValue("acceptable_obj_change_tol", 1e-4);
        //app->Options()->SetStringValue("hessian_approximation","limited-memory");
        //app->Options()->SetIntegerValue("print_level", 8);
        Ipopt::ApplicationReturnStatus status;
        status = app->Initialize();
        if (status != Ipopt::Solve_Succeeded) {
            std::cout << "# [error] Error during initialization!" << std::endl;
            return;
        }

        Ipopt::SmartPtr<Ipopt::TNLP> mynlp(
            new solver_use_ipopt(obj_, init_x));

        status = app->OptimizeTNLP(mynlp);

        if (status == Ipopt::Solve_Succeeded)
            std::cout << "# [info]  The problem solved!" << std::endl;
        else
            std::cout << "# [error] The problem failed!" << std::endl;

        return 0;
    }
private:
    third_solver_ipopt() {}
    func_ptr obj_;
};
}
}

#endif // THIRD_SOLVER_IPOPT_H
