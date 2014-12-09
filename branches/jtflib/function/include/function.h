#ifndef JTF_FUNCTION_H_
#define JTF_FUNCTION_H_

#include <set>
#include <vector>
#include <type_traits>
#include <memory>
#include <iostream>
#include <boost/shared_ptr.hpp>

namespace jtf { namespace function {

    //! @brief: R^n -> R

    class functionN1
    {
    public:
      virtual ~functionN1(){}

      virtual size_t dim(void) const = 0;

    public: // template functions
      //NOTE: add to v
      template <typename VAL_TYPE, typename INT_TYPE>
      int val(const VAL_TYPE *x, VAL_TYPE &v) {
        static_assert(std::is_same<double,VAL_TYPE>::value
                      ||std::is_same<float,VAL_TYPE>::value,
                      "bad value type");
        return val(reinterpret_cast<const void *>(x),
                   reinterpret_cast<void *>(v));
      }

      //NOTE: add to v
      template <typename VAL_TYPE, typename INT_TYPE>
      int gra(const VAL_TYPE *x, VAL_TYPE *g) {
        static_assert(std::is_same<VAL_TYPE,double>::value
                      ||std::is_same<VAL_TYPE,float>::value,
                      "bad value type");
        return gra(reinterpret_cast<const void*>(x),
                   reinterpret_cast<void*>(g));
      }


      template <typename VAL_TYPE, typename INT_TYPE>
      int gra(const VAL_TYPE *x, size_t &nnz, VAL_TYPE *g, INT_TYPE *idx)
      {return __LINE__;}


      //NOTE: add to hes
      //! @param h,ptr,idx == 0 means query nnz
      //! @param h == 0, ptr,idx != 0 means query patten
      //! @param h,ptr,idx != 0 means accumulate
      //! @param format for query pattern: 1 csc, 2 pair of row, col
      template <typename VAL_TYPE, typename INT_TYPE>
      int hes(const VAL_TYPE *x, size_t &nnz, size_t &format, VAL_TYPE *h,
              INT_TYPE *ptr, INT_TYPE *idx, VAL_TYPE alpha = 1){
        return hes(reinterpret_cast<const void*>(x),
                   nnz, format, reinterpret_cast<void*>(h),
                   reinterpret_cast<void*>(ptr), reinterpret_cast<void*>(idx), alpha);
      }


      template <typename VAL_TYPE, typename INT_TYPE>
      int hes_block(const VAL_TYPE *x, VAL_TYPE *h, VAL_TYPE alpha=1);//{ assert(0); return -1; }


      template <typename VAL_TYPE, typename INT_TYPE>
      bool is_valid(const VAL_TYPE *x) const {return true;}

    protected:
      //NOTE: add to v
      virtual int val(const void *x, void *v) = 0;

      //NOTE: add to v
      virtual int gra(const void *x, void *g) = 0;

      virtual int gra(const void *x, size_t &nnz, void *g, void *idx)
      {return __LINE__;}


      //NOTE: add to hes
      //! @param h,ptr,idx == 0 means query nnz
      //! @param h == 0, ptr,idx != 0 means query patten
      //! @param h,ptr,idx != 0 means accumulate
      //! @param format for query pattern: 1 csc, 2 pair of row, col
      virtual int hes(const void *x, size_t &nnz, size_t &format, void *h,
                      void *ptr, void *idx, double alpha = 1) = 0;

      virtual int hes_block(const void *x, void *h, double alpha=1) {return -1; }

      virtual bool is_valid(const void *x) const {return true;}
    };

    //! @brief: R^n -> R
    template <typename VAL_TYPE, typename INT_TYPE>
    class functionN1_t : public functionN1
    {
    public:
      typedef VAL_TYPE val_type;
      typedef INT_TYPE int_type;

      virtual ~functionN1_t(){}

      virtual size_t dim(void) const = 0;
      //NOTE: add to v
      virtual int val(const VAL_TYPE *x, VAL_TYPE &v) = 0;
      //NOTE: add to v
      virtual int gra(const VAL_TYPE *x, VAL_TYPE *g) = 0;
      virtual int gra(const VAL_TYPE *x, size_t &nnz, VAL_TYPE *g, INT_TYPE *idx) = 0;
      //NOTE: add to hes
      //! @param h,ptr,idx == 0 means query nnz
      //! @param h == 0, ptr,idx != 0 means query patten
      //! @param h,ptr,idx != 0 means accumulate
      //! @param format for query pattern: 1 csc, 2 pair of row, col
      virtual int hes(const VAL_TYPE *x, size_t &nnz, size_t &format, VAL_TYPE *h,
                      INT_TYPE *ptr, INT_TYPE *idx, VAL_TYPE alpha = 1) = 0;
      virtual int hes_block(const VAL_TYPE *x, VAL_TYPE *h, VAL_TYPE alpha=1) = 0 ;//{ assert(0); return -1; }
      virtual bool is_valid(const VAL_TYPE *x) const {return true;}

    protected:
      virtual int val(const void *x, void *v)  {
        return val(reinterpret_cast<const VAL_TYPE*>(x), *(reinterpret_cast<VAL_TYPE*>(v)));
      }

      //NOTE: add to v
      virtual int gra(const void *x, void *g){
        return gra(reinterpret_cast<const VAL_TYPE*>(x), reinterpret_cast<VAL_TYPE*>(g));
      }

      virtual int gra(const void *x, size_t &nnz, void *g, void *idx)
      {return gra(reinterpret_cast<const VAL_TYPE*>(x), nnz, reinterpret_cast<VAL_TYPE*>(g), reinterpret_cast<INT_TYPE*>(idx));}


      //NOTE: add to hes
      //! @param h,ptr,idx == 0 means query nnz
      //! @param h == 0, ptr,idx != 0 means query patten
      //! @param h,ptr,idx != 0 means accumulate
      //! @param format for query pattern: 1 csc, 2 pair of row, col
      virtual int hes(const void *x, size_t &nnz, size_t &format, void *h,
                      void *ptr, void *idx, double alpha = 1){
        return hes(reinterpret_cast<const VAL_TYPE*>(x), nnz, format, reinterpret_cast<VAL_TYPE*>(h),
                   reinterpret_cast<INT_TYPE*>(ptr), reinterpret_cast<INT_TYPE*>(idx), alpha);
      }

      virtual int hes_block(const void *x, void *h, double alpha=1) {
        return hes_block(reinterpret_cast<const VAL_TYPE*>(x), reinterpret_cast<VAL_TYPE*>(h), alpha);}

      virtual bool is_valid(const void *x) const {
        return is_valid(reinterpret_cast<const VAL_TYPE*>(x));}
    };

    enum POINT_TYPE{RAW, SMART_BOOST,SMART_BOOST_CONS, SMART_STD,SMART_STD_CONS};

    template <typename VAL_TYPE, typename INT_TYPE, POINT_TYPE PT>
    class point_type_traits;

    template <typename VAL_TYPE, typename INT_TYPE>
    class point_type_traits<VAL_TYPE, INT_TYPE, RAW>{
    public:
      typedef std::vector<functionN1_t<VAL_TYPE,INT_TYPE> const*> container;      
    };

    template <typename VAL_TYPE, typename INT_TYPE>
    class point_type_traits<VAL_TYPE, INT_TYPE, SMART_STD_CONS>{
    public:
      typedef std::vector<std::shared_ptr<const functionN1_t<VAL_TYPE,INT_TYPE> > > container;
    };

    template <typename VAL_TYPE, typename INT_TYPE>
    class point_type_traits<VAL_TYPE, INT_TYPE, SMART_STD>{
    public:
      typedef std::vector<std::shared_ptr<functionN1_t<VAL_TYPE,INT_TYPE> > > container;
    };

    template <typename VAL_TYPE, typename INT_TYPE>
    class point_type_traits<VAL_TYPE, INT_TYPE, SMART_BOOST_CONS>{
    public:
      typedef std::vector<boost::shared_ptr<const functionN1_t<VAL_TYPE,INT_TYPE> > > container;
    };

    template <typename VAL_TYPE, typename INT_TYPE>
    class point_type_traits<VAL_TYPE, INT_TYPE, SMART_BOOST>{
    public:
      typedef std::vector<boost::shared_ptr<functionN1_t<VAL_TYPE,INT_TYPE> > > container;
    };

    template <typename VAL_TYPE, typename INT_TYPE, POINT_TYPE PT = SMART_STD>
    class sum_function : public functionN1_t<VAL_TYPE, INT_TYPE>
    {
    public:
      typedef typename point_type_traits<VAL_TYPE, INT_TYPE, PT>::container container;
      sum_function(const container &children):children_(children), nnz_(-1)
      {
        const size_t fn = children_.size();
        for(size_t i = 1; i < fn; ++i) {
            if(children_[i]->dim() != children_[0]->dim()) {
                std::cerr << "incompatible functions." << children_[0]->dim()
                          << " " << children_[i]->dim() << std::endl;
              }
          }
      }
      virtual size_t dim(void) const {return children_[0]->dim();}
      virtual int val(const VAL_TYPE *x, VAL_TYPE &v){
        const size_t fn = children_.size();
        for(const auto & func_ptr : children_) {
            func_ptr->val(x,v);
          }

        return 0;
      }
      virtual int gra(const VAL_TYPE *x, VAL_TYPE *g){
        for(const auto & func_ptr : children_) func_ptr->gra(x,g);
        return 0;
      }
      virtual int gra(const VAL_TYPE *x, size_t &nnz, VAL_TYPE *g, INT_TYPE *idx){

        std::set<INT_TYPE> idx_;
        std::vector<INT_TYPE> idx_vec;
        std::vector<VAL_TYPE> g_vec;
        std::vector<VAL_TYPE> full_g(dim());
        fill(full_g.begin(), full_g.end(), 0);
        for(const auto & func_ptr: children_){
            func_ptr->gra(x, nnz, 0, 0);
            idx_vec.resize(nnz);
            g_vec.resize(nnz);

            func_ptr->gra(x,nnz, &g_vec[0], &idx_vec[0]);
            for(size_t i = 0; i < idx_vec.size(); ++i)
              full_g[idx_vec[i]] += g_vec[i];

            idx_.insert(idx_vec.begin(), idx_vec.end());
          }

        if(g == 0 && idx == 0)
          nnz = idx_.size();
        if(g != 0 && idx != 0){
            std::copy(idx_.begin(), idx_.end(), idx);
            size_t idx_i = 0;
            for(const auto & one_idx: idx_)
              g[idx_i++] = full_g[one_idx];
          }

        return 0;
      }

      virtual int hes(const VAL_TYPE *x, size_t &nnz, size_t &format, VAL_TYPE *h, INT_TYPE *ptr, INT_TYPE *idx, double alpha = 1)
      {
        format = 1;
        if(h == 0 && ptr == 0 && idx == 0) {// query nnz
            pattern_.resize(dim());
            std::pair<std::vector<INT_TYPE>, std::vector<INT_TYPE> > ptr_idx;
            for(size_t fi = 0; fi < children_.size(); ++fi) {
                size_t nnz0, format = -1;
                if(children_[fi]->hes(x, nnz0, format, 0, 0, 0))
                  return __LINE__;
                if(format == 1) { // csc
                    ptr_idx.first.clear();
                    ptr_idx.first.resize(dim()+1);
                    ptr_idx.first[0] = 0;
                    ptr_idx.second.clear();
                    ptr_idx.second.resize(nnz0);
                    if(children_[fi]->hes(x, nnz0, format, 0, &ptr_idx.first[0], &ptr_idx.second[0]))
                      return __LINE__;
                    for(size_t ci = 0; ci < dim(); ++ci) {
                        for(size_t nzi = ptr_idx.first[ci]; nzi < ptr_idx.first[ci+1]; ++nzi) {
                            pattern_[ci].insert(ptr_idx.second[nzi]);
                          }
                      }
                  }
                else if(format == 2) {// pair
                    ptr_idx.first.resize(nnz0);
                    ptr_idx.second.resize(nnz0);
                    if(children_[fi]->hes(x, nnz0, format, 0, &ptr_idx.first[0], &ptr_idx.second[0]))
                      return __LINE__;
                    for(size_t nzi = 0; nzi < nnz0; ++nzi) {
                        pattern_[ptr_idx.first[nzi]].insert(ptr_idx.second[nzi]);
                      }
                  }
              }
            nnz = 0;
            for(size_t xi = 0; xi < dim(); ++xi)
              nnz += pattern_[xi].size();
            nnz_ = nnz;
            return 0;
          }
        if(h == 0 && ptr != 0 && idx != 0) {// query patten
            if(nnz < nnz_) {
                std::cerr << "incorrect input at query pattern: " << nnz << " " << nnz_;
                return __LINE__;
              }
            for(size_t xi = 0; xi < dim(); ++xi) {
                ptr[xi+1] = ptr[xi] + pattern_[xi].size();
                size_t nzi = ptr[xi];
                for(typename std::set<INT_TYPE>::const_iterator i = pattern_[xi].begin();
                    i != pattern_[xi].end(); ++i, ++nzi) {
                    idx[nzi] = *i;
                  }
              }
            std::vector<std::set<INT_TYPE> > tmp;
            std::swap(pattern_, tmp);
            return 0;
          }
        if(h != 0 && ptr != 0 && idx != 0) {// accumulate
            if(nnz < nnz_ && nnz_ != -1) { // when nnz_ == -1, client know the pattern already
                std::cerr << "incorrect input at accumulate: " << nnz << " " << nnz_;
                return __LINE__;
              }
            size_t format = -1;
            for(size_t fi = 0; fi < children_.size(); ++fi)
              if(children_[fi]->hes(x, nnz, format, h, ptr, idx, alpha))
                return __LINE__;
            return 0;
          }
        return __LINE__;
      }
      virtual int hes_block(const VAL_TYPE *x, VAL_TYPE *h, VAL_TYPE alpha = 1) {return -1;}
      virtual bool is_valid(const VAL_TYPE *x) const{
        for(size_t fi = 0; fi < children_.size(); ++fi)
          if(!children_[fi]->is_valid(x))
            return false;
        return true;
      }
    protected:
      const container children_;
      std::vector<std::set<INT_TYPE> > pattern_;
      size_t nnz_;
    };
  }}

#endif
