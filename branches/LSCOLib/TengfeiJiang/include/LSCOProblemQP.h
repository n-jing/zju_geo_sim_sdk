#ifndef _LSCO_PROBLEM_QP_H_
#define _LSCO_PROBLEM_QP_H_

#include "kernel_common_file_pzr.h"
#include "direct_solver_kernel_common_file_pzr.h"
#include "objective_function.h"
#include <set>

namespace zjucad
{
namespace LSCO
{

using namespace COMMON;
typedef kernel_traits<FixedSparseMatrix<scalarD,Kernel<scalarD> > > PROBLEM_KERNEL;

class ProblemQP : public objective_function<PROBLEM_KERNEL>
{
public:
    ProblemQP(sizeType nr,sizeType nrE,sizeType nrI,bool diagonal,bool feasible=true)
        :_nr(nr),_nrE(nrE),_nrI(nrI),_diagonal(diagonal) {
        _H.resize(nr,nr);
        _g.resize(nr);
        for(sizeType i=0; i<nr; i++)
            for(sizeType j=0; j<nr; j++) {
                if(j > i) {
                    value_type v=rand()*2.0f/(value_type)RAND_MAX-1.0f;
                    _H.addToElement(j,i,v);
                    _H.addToElement(i,j,v);
                }
                _g[i]=rand()*2.0f/(value_type)RAND_MAX-1.0f;
            }
        for(sizeType i=0; i<nr; i++) {
            value_type dv=20.0f;
            for(sizeType j=0; j<nr; j++)
                dv+=std::abs(_H(i,j));
            _H.addToElement(i,i,dv);
        }

		_IA.resize(nrE+nrI,nr);
		_bd.resize(nrE+nrI);

        _A.resize(nrE,nr);
        _b.resize(nrE);
        std::set<sizeType> visited;
        for(sizeType i=0; i<nrE; i++) {
            sizeType a=rand()%nr;
            while(visited.find(a) != visited.end())a=rand()%nr;
            visited.insert(a);

            sizeType b=rand()%nr;
            while(visited.find(b) != visited.end())b=rand()%nr;
            visited.insert(b);

            value_type ca=rand()*2.0f/(value_type)RAND_MAX-1.0f;
            value_type cb=rand()*2.0f/(value_type)RAND_MAX-1.0f;
            value_type cv=rand()*2.0f/(value_type)RAND_MAX-1.0f;
            _A.addToElement(i,a,ca);
            _A.addToElement(i,b,cb);
			_IA.addToElement(i,a,ca);
            _IA.addToElement(i,b,cb);
            _b[i]=cv;
			_bd[i]=cv;
        }

        visited.clear();
        _I.resize(nrI,nr);
        _d.resize(nrI);
        for(sizeType i=0; i<nrI; i++) {
            sizeType a=rand()%nr;
            while(visited.find(a) != visited.end())a=rand()%nr;
            visited.insert(a);

            sizeType b=rand()%nr;
            while(visited.find(b) != visited.end())b=rand()%nr;
            visited.insert(b);

            value_type ca=rand()*2.0f/(value_type)RAND_MAX-1.0f;
            value_type cb=rand()*2.0f/(value_type)RAND_MAX-1.0f;
            value_type cv=rand()*2.0f/(value_type)RAND_MAX-1.0f;
            _I.addToElement(i,a,ca);
            _I.addToElement(i,b,cb);
			_IA.addToElement(nrE+i,a,ca);
            _IA.addToElement(nrE+i,b,cb);
            _d[i]=cv;
			_bd[nrE+i]=cv;
        }

		if(feasible) {
            boost::property_tree::ptree pt;
            direct_solver_kernel<PROBLEM_KERNEL::sparse_matrix_type> sol(pt);
            sol.build_AAT(_IA,false,false);

            vector_type rhs;
            rhs.resize(_nr);
            rhs.setZero();
            _IA.multiplyTransposeSubtract(_bd,rhs);
            sol.solve(rhs,_x);
        }
    }
    virtual size_t nr_inputs() const {
        return _nr;
    }
    virtual size_t nr_equals() const {
        return _nrE;
    }
    virtual size_t nr_inequals() const {
        return _nrI;
    }
    virtual value_type eval(const vector_type& x,bool same_x) {
        if(same_x)
            assert((_old_x-x.block(0,0,_nr,1)).norm() == 0.0f);
        else _old_x=x.block(0,0,_nr,1);

        vector_type Hx=_g;
        _H.multiply(x,Hx);
        return PROBLEM_KERNEL::dot(x,Hx,nr_inputs())*0.5f+PROBLEM_KERNEL::dot(x,_g,nr_inputs())+10.0f;
    }
    virtual void eval_func_gradient(const vector_type& x,vector_type& gradient,bool same_x) {
        if(same_x)
            assert((_old_x-x.block(0,0,_nr,1)).norm() == 0.0f);
        else _old_x=x.block(0,0,_nr,1);

        _H.multiply(x,gradient);
        PROBLEM_KERNEL::axpy(1.0f,_g,gradient,nr_inputs());
    }
    virtual void eval_lagrangian_hessian(const vector_type& x,const vector_type& lambda,value_type coef,PROBLEM_KERNEL::matrix_builder& hessian,bool same_x) {
        if(same_x)
            assert((_old_x-x.block(0,0,_nr,1)).norm() == 0.0f);
        else _old_x=x.block(0,0,_nr,1);

        if(_diagonal) {
            hessian.clear();
            for(sizeType r=0; r<_H.rows(); r++)
                hessian.push_back(Eigen::Triplet<value_type,sizeType>(r,r,1.0f));
        } else {
            hessian.clear();
            for(sizeType r=0; r<_H.rows(); r++)
                for(ConstSMIterator<value_type> beg=_H.begin(r),end=_H.end(r); beg!=end; ++beg)
                    hessian.push_back(Eigen::Triplet<value_type,sizeType>(r,beg.col(),*beg));
        }
    }
    virtual void eval_constraint_value(const vector_type& x,vector_type& c,bool same_x) {
        if(same_x)
            assert((_old_x-x.block(0,0,_nr,1)).norm() == 0.0f);
        else _old_x=x.block(0,0,_nr,1);

        vector_type EV;
        vector_type IV;
        EV.resize(_nrE);
        IV.resize(_nrI);
        _A.multiply(x,EV);
        _I.multiply(x,IV);
        EV+=_b;
        IV+=_d;
        c.block(0,0,_nrE,1)=EV;
        c.block(_nrE,0,_nrI,1)=IV;
    }
    virtual void eval_constraint_gradient(const vector_type& x,PROBLEM_KERNEL::matrix_builder& c_gradient,bool same_x) {
        if(same_x)
            assert((_old_x-x.block(0,0,_nr,1)).norm() == 0.0f);
        else _old_x=x.block(0,0,_nr,1);

        c_gradient.clear();
        for(sizeType r=0; r<_A.rows(); r++)
            for(ConstSMIterator<value_type> beg=_A.begin(r),end=_A.end(r); beg!=end; ++beg)
                c_gradient.push_back(Eigen::Triplet<value_type,sizeType>(r,beg.col(),*beg));
        for(sizeType r=0; r<_I.rows(); r++)
            for(ConstSMIterator<value_type> beg=_I.begin(r),end=_I.end(r); beg!=end; ++beg)
                c_gradient.push_back(Eigen::Triplet<value_type,sizeType>(r+_nrE,beg.col(),*beg));
    }
    vector_type _x;
protected:
    sizeType _nr,_nrE,_nrI;
    sparse_matrix_type _H,_A,_I,_IA;
    vector_type _g,_b,_d,_bd,_old_x;
    bool _diagonal;
};

}
}

#endif