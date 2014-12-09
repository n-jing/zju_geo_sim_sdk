#ifndef JTF_PRECONDITION_H
#define JTF_PRECONDITION_H

#include <zjucad/matrix/matrix.h>
#include <hjlib/sparse/sparse.h>

namespace jtf{
  enum PRE_TYPE{LEFT,RIGHT};


  class precondition{
  public:
    virtual ~precondition(){}
  };

  class jacobi_precondition: public precondition
  {
  public:
    jacobi_precondition(const zjucad::matrix::matrix<double> & A,
                        const zjucad::matrix::matrix<double> & b)
      :A_(A), b_(b), weight_(1.0){init();}

    virtual ~jacobi_precondition(){}
  public:
    virtual void get_r(const double *x, double *r) const;

    virtual void get_b(double *b) const;

    virtual void get_Ax(const double *x, double *Ax)const;
  private:
    void init();
    const zjucad::matrix::matrix<double> & A_;
    const zjucad::matrix::matrix<double> & b_;
    zjucad::matrix::matrix<double> M_;
    const double weight_;
  };

  class jacobi_precondition_sparse: public precondition
  {
  public:
    jacobi_precondition_sparse(const hj::sparse::csc<double,int32_t> & A,
                               const zjucad::matrix::matrix<double> & b)
      :A_(A), b_(b), weight_(1.0){init();}

    virtual ~jacobi_precondition_sparse(){}
  public:
    virtual void get_r(const double *x, double *r) const;

    virtual void get_b(double *b) const;

    virtual void get_Ax(const double *x, double *Ax)const;
  private:
    void init();
    const hj::sparse::csc<double,int32_t> & A_;
    const zjucad::matrix::matrix<double> & b_;
    zjucad::matrix::matrix<double> M_;
    const double weight_;
  };

}

#endif
