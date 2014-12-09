#ifndef OPERATION_H
#define OPERATION_H

#include "function.h"

namespace jtf{
  namespace  function {

    template <template <typename VAL_TYPE> class OP, typename VAL_TYPE,
              typename INT_TYPE, template <typename FUNC> class PTR>
    class scalar_op : public jtf::function::functionN1_t<VAL_TYPE, INT_TYPE>
    {
    public:
      typedef VAL_TYPE value_type;
      typedef INT_TYPE int_type;

      scalar_op(const functionN1_t<value_type,int_type> &src, value_type scalar)
        :src_(src), scalar_(scalar) {
      }
      scalar_op(PTR<const functionN1_t<value_type,int_type>> own, value_type scalar)
        :src_(*own), own_(own), scalar_(scalar) {
      }
      virtual size_t dim(void) const {return src_.dim();}
      virtual int val(const value_type *x, value_type &f) const {
        value_type f_temp = 0;
        if(src_.val(x, f_temp))
          return __LINE__;
        f += OP<VAL_TYPE>()(f_temp, scalar_);
        return 0;
      }
      virtual int gra(const value_type *x, value_type *g){
        static std::vector<value_type> g_empty;
        if(g_empty.size() != dim()) g_empty.resize(dim(),0);
        std::fill(g_empty.begin(), g_empty.end(),0);

        if(src_.gra(x,&g_empty[0])) return __LINE__;
        for(size_t i = 0; i < dim(); ++i)
          g[i] += OP<VAL_TYPE>()(g_empty[i], scalar_);
        return 0;
      }
      virtual int gra(const value_type *x, size_t &nnz, value_type * g, int32_t *idx){

        if(g == 0 && idx == 0){
            if(src_.gra(x, nnz, g, idx)) return __LINE__;
          }else{
            static std::vector<value_type> g_empty;
            if(g_empty.size() != nnz) g_empty.resize(nnz, 0);
            std::fill(g_empty.begin(), g_empty.end(),0);
            if(src_.gra(x, nnz, &g_empty[0], idx)) return __LINE__;

            for(size_t i = 0; i < nnz; ++i) g[i] += OP<VAL_TYPE>()(g_empty[i],scalar_);
          }
        return 0;
      }
      virtual int hes(const value_type *x, size_t &nnz, size_t &format, value_type *h,
                      int_type *ptr, int_type *idx, double alpha = 1)
      {
        if(h == 0) {
            return src_.hes(x, nnz, format, h, ptr, idx, alpha);
          }
        static std::vector<value_type> h_empty;
        if(h_empty.size() != nnz) h_empty.resize(nnz, 0);
        std::fill(h_empty.begin(), h_empty.end(),0);
         src_.hes(x, nnz, format, &h_empty[0], ptr, idx, alpha);
        for(size_t i = 0; i < nnz; ++i) h[i] += OP<value_type>()(h_empty[i], scalar_);
        return 0;
      }
      virtual int hes_block(const value_type *x, value_type *h, double alpha)
      {
        throw "# [error] I do not know how to calculate hes_block.";
      }
    protected:
      value_type scalar_;
      const jtf::function::functionN1_t<VAL_TYPE,INT_TYPE> &src_;
      PTR<const jtf::function::functionN1_t<VAL_TYPE,INT_TYPE> > own_;
    };

    template <typename FUNC>
    class nil_ptr {};

#define JTF_SCALAR_OP_REF(OP, NAME) \
  template <typename VAL_TYPE, typename INT_TYPE> \
  jtf::function::functionN1_t<VAL_TYPE, INT_TYPE> \
  *operator OP (const jtf::function::functionN1_t<VAL_TYPE, INT_TYPE> &func, VAL_TYPE scalar) { \
  return new scalar_op<NAME, VAL_TYPE, INT_TYPE, nil_ptr>(func, scalar);		\
  }

    JTF_SCALAR_OP_REF(+, std::plus)
    JTF_SCALAR_OP_REF(-, std::minus)
    JTF_SCALAR_OP_REF(*, std::multiplies)
    JTF_SCALAR_OP_REF(/, std::divides)


#define JTF_SCALAR_OP_OWN(OP, NAME) \
  template <typename VAL_TYPE, typename INT_TYPE,template <typename FUNC> class PTR> \
  jtf::function::functionN1_t<VAL_TYPE, INT_TYPE> \
  *operator OP (const PTR<const functionN1_t<VAL_TYPE, INT_TYPE> > func, VAL_TYPE scalar) { \
  return new scalar_op<NAME, VAL_TYPE, INT_TYPE, PTR>(func, scalar);		\
  }

    JTF_SCALAR_OP_OWN(+, std::plus)
    JTF_SCALAR_OP_OWN(-, std::minus)
    JTF_SCALAR_OP_OWN(*, std::multiplies)
    JTF_SCALAR_OP_OWN(/, std::divides)
  }
}
#endif // OPERATION_H
