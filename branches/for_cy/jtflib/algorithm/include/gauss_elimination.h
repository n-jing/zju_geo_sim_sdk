#ifndef GAUSS_ELIMINATION_H
#define GAUSS_ELIMINATION_H

#include <map>
#include <set>
#include <vector>
#include <list>
#include <deque>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/unordered_map.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>
#include <zjucad/matrix/matrix.h>

#include "equation.h"

namespace jtf {
  namespace algorithm{

    //! @brief this class only handle Ai+Bi=Ci
    template <typename T>
    class gauss_eliminator{
      BOOST_MPL_ASSERT_MSG((boost::is_same<T,double>::value ) ||
                           (boost::is_same<T,float>::value ),
                           NON_FLOAT_TYPES_ARE_NOT_ALLOWED, (void));

    public:
      typedef typename std::list<equation<T> >::iterator equation_ptr;
      typedef typename std::list<equation<T> >::const_iterator const_equation_ptr;
      typedef typename std::map<size_t, std::list<equation_ptr> >::iterator prime_eq_ptr;

      /**
 * @brief construct gauss_eliminator class
 *
 * @param nodes input nodes
 * @param node_flag input node_flag which will be tagged as true if the
 *        corresponding node is known
 */
      explicit gauss_eliminator(
          std::vector<T> & nodes,
          boost::dynamic_bitset<> & node_flag,
          bool do_nothing = false)
      //std::vector<bool> & node_flag)
        :nodes_(nodes), node_flag_(node_flag), do_nothing_(do_nothing){
        idx2equation_.resize(nodes_.size());
      }

      const_equation_ptr begin()const{return es.begin();}
      const_equation_ptr end()const{return es.end();}

      int clear(){
        es.clear();
        idx2equation_.clear();
        idx2equation_.resize(nodes_.size());
        prime_idx2equation_.clear();
      }

      equation_ptr begin(){return es.begin();}
      equation_ptr end(){return es.end();}

      int add_equation(const equation<T> & e);

      size_t get_eqn_number() const {return es.size();}
      // convert to A*X = B
      int convert_to_matrix(zjucad::matrix::matrix<T> & A,
                            zjucad::matrix::matrix<T> & B)const;

      const std::list<equation<T> > & get_equations()const {return es;}

      int to_octave_sparse_matrix_triples(const std::string &file_name )const;

      int save_equations(std::ostream & ofs)const;
      int load_equations(std::istream & ifs);

      int save_equations_mem(std::vector<T> & node_vec,
                             boost::dynamic_bitset<> & node_flag,
                             std::vector<equation<T> > & eq_vec)const;

      int load_equations_mem(const std::vector<T> & node_vec,
                             const boost::dynamic_bitset<> & node_flag,
                             const std::vector<equation<T> > & eq_vec);

      bool check_gauss_eliminator()
      {
        std::set<size_t> prim_index;
        for(const auto & one_eqn : es){
            size_t prim_index_i = one_eqn.get_prime_idx();
            const size_t number = prim_index.size();
            prim_index.insert(prim_index_i);
            if(prim_index.size() == number){
                std::cerr << "# [error] two equation with the same prim index." << std::endl;
                return false;
              }
            auto it = one_eqn.begin();
            if(it->index != prim_index_i){
                std::cerr << "# [error] fisrst term is not prim term." << std::endl;
                return false;
              }
            ++it;
            for(; it != one_eqn.end(); ++it){
                if(std::find(prim_index.begin(), prim_index.end(),
                             it->index) != prim_index.end()){
                    std::cerr << "# [error] equation has other prim index."
                              << it->index << std::endl;
                    return false;
                  }
              }
          }
        return true;
      }

    private:
      int remove_idx2eqn(const equation_ptr eqn_it,
                         const size_t exclue_idx);

      int add_idx2eqn(const equation_ptr eqn_it,
                      const size_t exclue_idx);

      int remove_all_other_idx2eqn(const equation_ptr eqn_it,
                                   const size_t exclue_idx);

      int add_all_other_idx2eqn(const equation_ptr eqn_it,
                                const size_t exclue_idx);

      T get_coefficient(const equation_ptr eqn_it,
                        const size_t idx);

    private:
      std::vector<T> & nodes_;
      boost::dynamic_bitset<> & node_flag_;
      std::list<equation<T> > es;
      bool do_nothing_;
      std::vector<std::list<equation_ptr> > idx2equation_;

      // this map store the smallest expression
      std::map<size_t, std::list<equation_ptr> > prime_idx2equation_;
    private:
      gauss_eliminator(const gauss_eliminator<T>&){}
      gauss_eliminator<T>& operator=(const gauss_eliminator<T>&){}
      int gather_variant_index_of_equation(const equation<T> & eq,
                                           std::vector<size_t> & gather)const;
    };

    template <typename T>
    int gauss_eliminator<T>::add_equation(const equation<T> & e)
    {
      es.push_back(e);
      auto & e_back = es.back();
      for(typename equation<T>::eq_iterator eit = e_back.begin(); eit != e_back.end(); ){
          if(node_flag_[eit->index]){
              e_back.value() -= nodes_[eit->index] * eit->coefficient;
              e_back.e_vec_.erase(eit++);
            }else
            ++eit;
        }
      e_back.standardization();

      if(e_back.state() == 0){// this equation is cleared
          es.pop_back();
          return 0;
        }else if(e_back.state() == -1){
          std::cerr << "# [error] strange conflict equation: " << std::endl;
          std::cerr << e;
          es.pop_back();
          return __LINE__;
        }

      equation<T> e_temp;
      e_temp.value() = e_back.value();
      for(const auto & one_exp : e_back){
          const auto & it = prime_idx2equation_.find(one_exp.index);
          if(it != prime_idx2equation_.end()){
              const auto &temp = *(it->second.front());
              for(const auto & one_exp_temp : temp){
                  if(one_exp_temp.index == one_exp.index) continue;
                  e_temp.add_expression(
                        make_expression(
                          one_exp_temp.index,
                          -1*one_exp.coefficient * one_exp_temp.coefficient));
                }
              e_temp.value() += -1 * one_exp.coefficient * temp.value();
            }else{
              e_temp.add_expression(one_exp);
            }
        }
      e_temp.standardization();
      e_back = e_temp;

      if(e_back.state() == 0){// this equation is cleared
          es.pop_back();
          return 0;
        }else if(e_back.state() == -1){
          std::cerr << "# [error] strange conflict equation: " << std::endl;
          std::cerr << e;
          es.pop_back();
          return __LINE__;
        }

      // this equation can be inserted into es, these stored equations may be eliminated by it.
      const size_t prime_idx = e_back.get_prime_idx();
      auto it_this_eqn = es.end();
      --it_this_eqn;
      prime_idx2equation_[prime_idx].push_back(it_this_eqn);

      if(e_back.state() == 1){
          node_flag_[prime_idx] = true;
          nodes_[prime_idx] = e_back.value();
        }

      auto & idx2other_eqn = idx2equation_[prime_idx];
      std::set<size_t> idx_before,idx_after;
      std::vector<size_t> temp_container;
      bool is_erased = false;
      for(auto it  = idx2other_eqn.begin(); it != idx2other_eqn.end();++it){
          is_erased = false;
          double w = get_coefficient(*it, prime_idx);
          idx_before.clear();
          for(const auto & one_exp : *(*it))
            idx_before.insert(one_exp.index);

          if(fabs(w-1) <1e-6 ){
              *(*it) -= e_back;
            }else{
              auto temp = e_back;
              temp *= w;
              *(*it) -= temp;
            }
          if((*it)->state() == 0) {
              es.erase(*it);
              is_erased = true;
            }else if((*it)->state() == 1){
              const size_t prime_idx_1 = (*it)->get_prime_idx();
              node_flag_[prime_idx_1] = true;
              nodes_[prime_idx_1] = (*it)->value();
            }
          {
            idx_after.clear();
            if(is_erased == false){
                for(const auto & one_exp : *(*it))
                  idx_after.insert(one_exp.index);
              }
            temp_container.resize(std::max(idx_before.size(), idx_after.size()));
            auto cit_need_to_remove_end =
                std::set_difference(idx_before.begin(),idx_before.end(),
                                    idx_after.begin(),idx_after.end(),
                                    temp_container.begin());
            for(auto it_idx = temp_container.begin();
                it_idx != cit_need_to_remove_end; ++it_idx){
                if(*it_idx == prime_idx) continue;
                remove_idx2eqn(*it, *it_idx);
              }

            if(is_erased == false){
                auto cit_need_to_add_end =
                    std::set_difference(idx_after.begin(),idx_after.end(),
                                        idx_before.begin(),idx_before.end(),
                                        temp_container.begin());
                for(auto it_idx = temp_container.begin(); it_idx != cit_need_to_add_end; ++it_idx){
                    if(*it_idx == prime_idx) continue;
                    add_idx2eqn(*it, *it_idx);
                  }
              }
          }
        }

      idx2other_eqn.clear();

      add_all_other_idx2eqn(it_this_eqn, prime_idx);
    }

    template <typename T>
    T gauss_eliminator<T>::get_coefficient(const equation_ptr eqn_it,
                                           const size_t idx)
    {
      for(const auto & one_exp : *eqn_it){
          if(one_exp.index == idx)
            return one_exp.coefficient;
        }
      std::cerr << "# [error] strange, can not find index " << idx << " in eqn." << std::endl;
      return 0;
    }

    template <typename T>
    int gauss_eliminator<T>::remove_all_other_idx2eqn(const equation_ptr eqn_it,
                                                      const size_t exclue_idx)
    {
      const size_t prime_idx = eqn_it->get_prime_idx();
      for(const auto & one_exp : *eqn_it){
          if(exclue_idx == one_exp.index
             || prime_idx == one_exp.index) continue;
          auto & list_of_eqn = idx2equation_[one_exp.index];
          auto it = std::find(list_of_eqn.begin(), list_of_eqn.end(), eqn_it);
          assert(it != list_of_eqn.end());
          list_of_eqn.erase(it);
        }
      return 0;
    }

    template <typename T>
    int gauss_eliminator<T>::remove_idx2eqn(const equation_ptr eqn_it,
                                            const size_t exclue_idx)
    {
      auto & list_of_eqn = idx2equation_[exclue_idx];
      auto it = std::find(list_of_eqn.begin(), list_of_eqn.end(), eqn_it);
      if(it != list_of_eqn.end())
        list_of_eqn.erase(it);

      return 0;
    }

    template <typename T>
    int gauss_eliminator<T>::add_idx2eqn(const equation_ptr eqn_it,
                                         const size_t exclue_idx)
    {
      auto & list_of_eqn = idx2equation_[exclue_idx];
      // speed up !
      //      auto it = std::find(list_of_eqn.begin(), list_of_eqn.end(), eqn_it);
      //      assert(it == list_of_eqn.end());
      list_of_eqn.push_back(eqn_it);

      return 0;
    }

    template <typename T>
    int gauss_eliminator<T>::add_all_other_idx2eqn(const equation_ptr eqn_it,
                                                   const size_t exclue_idx)
    {
      const size_t prime_idx = eqn_it->get_prime_idx();
      for(const auto & one_exp : *eqn_it){
          if(exclue_idx == one_exp.index ||
             prime_idx == one_exp.index) continue;
          auto & list_of_eqn = idx2equation_[one_exp.index];
          auto it = std::find(list_of_eqn.begin(), list_of_eqn.end(), eqn_it);
          assert(it == list_of_eqn.end());
          list_of_eqn.push_back(eqn_it);
        }
      return 0;
    }

    template <typename T>
    int gauss_eliminator<T>::to_octave_sparse_matrix_triples(
        const std::string & filename)const
    {
      std::ofstream ofs(filename.c_str());
      if(ofs.fail()){
          std::cerr << " can not open octave file." << std::endl;
          return __LINE__;
        }

      ofs << "# Created by jtf::gauss_elimination" << std::endl;
      ofs << "# name: A" << std::endl;
      ofs << "# type: sparse matrix" << std::endl;

      std::vector<std::tuple<size_t,size_t,T> > triples;
      size_t rows_num = 0;
      size_t nnz = 0;

      for(typename std::list<equation<T> >::const_iterator lcit = es.begin();
          lcit != es.end(); ++lcit){
          if((*lcit).state() == 0) continue;
          if((*lcit).state() == -1) continue;

          const equation<T> & eq = *lcit;
          const size_t triples_num = triples.size();
          for(typename equation<T>::eq_const_iterator eqcit = eq.begin();
              eqcit != eq.end(); ++eqcit){
              const expression<T> &exp = *eqcit;
              if(fabs(exp.coefficient) > 1e-6){
                  triples.push_back(std::make_tuple(
                                      rows_num, exp.index, exp.coefficient));
                  ++nnz;
                }
            }
          if(triples.size() != triples_num)
            ++rows_num;
        }

      ofs << "# nnz: " << nnz << std::endl;
      ofs << "# rows: " << rows_num << std::endl;
      ofs << "# columns: " << nodes_.size() << std::endl;

      for(size_t ni = 0; ni < triples.size(); ++ni){
          const std::tuple<size_t,size_t,T> & trip = triples[ni];
          ofs << std::get<0>(trip) << " " << std::get<1>(trip)
              << " " << std::get<2>(trip) << std::endl;
        }
      return 0;
    }


    template <typename T>
    int gauss_eliminator<T>::save_equations(std::ostream & ofs)const
    {
      //std::ofstream ofs(file_name.c_str(), std::ios::binary);
      if(ofs.fail()){
          std::cerr << "# [error] can not open equation file "  << std::endl;
          return __LINE__;
        }

      const size_t node_num = nodes_.size();
      ofs.write((char*)&node_num, sizeof(size_t));
      ofs.write((char *)&nodes_[0], nodes_.size() * sizeof(size_t));

      const size_t block_num = node_flag_.num_blocks();
      ofs.write((char*)&block_num, sizeof(size_t));
      std::vector<boost::dynamic_bitset<>::block_type> block_vec(block_num);
      boost::to_block_range(node_flag_, block_vec.begin());
      ofs.write((char*)(&block_vec[0]), block_num *
          sizeof(boost::dynamic_bitset<>::block_type));

      std::vector<const_equation_ptr> left_equations;
      left_equations.reserve(es.size());
      for(const_equation_ptr cep = es.begin(); cep != es.end(); ++cep){
          if(cep->state() == 0) continue;
          if(cep->state() == 1){
              // here takes an assumption, this equations is resolved
              continue;
            }
          if(cep->state() == -1){
              std::cerr << "# [error] this equation is conflict " << std::endl;
              std::cerr << *cep << std::endl;
              continue;
            }
          left_equations.push_back(cep);
        }
      const size_t eq_num = left_equations.size();
      ofs.write((char*)&eq_num, sizeof(size_t));
      for(size_t ei = 0; ei < left_equations.size(); ++ei){
          const equation<T> & eq = *(left_equations[ei]);
          eq.write(ofs);
        }

      return 0;
    }

    template <typename T>
    int gauss_eliminator<T>::load_equations(std::istream & ifs)
    {
      //std::ifstream ifs(file_name.c_str(), std::ios::binary);
      if(ifs.fail()){
          std::cerr << "# [error] can not open equation file." << std::endl;
          return __LINE__;
        }

      size_t node_num = 0;
      ifs.read((char*)&node_num, sizeof(size_t));
      nodes_.resize(node_num);
      ifs.read((char *)&nodes_[0], node_num * sizeof(size_t));

      size_t block_num = 0;
      ifs.read((char*)&block_num, sizeof(size_t));
      std::vector<boost::dynamic_bitset<>::block_type> block_vec(block_num);
      ifs.read((char*)&block_vec[0], block_num *
          sizeof(boost::dynamic_bitset<>::block_type));
      node_flag_.resize(node_num);
      boost::from_block_range(block_vec.begin(), block_vec.end(), node_flag_);

      size_t equation_num = 0;
      ifs.read((char*)&equation_num, sizeof(size_t));

      clear();

      for(size_t eqi = 0; eqi < equation_num; ++eqi){
          equation<T> eq;
          eq.read(ifs);
          es.push_back(eq);

          equation_ptr ep = es.end();
          --ep;
          prime_idx2equation_[eq.e_vec_.front().index].push_back(ep);

          for(typename equation<T>::eq_const_iterator eqcit = eq.begin();
              eqcit != eq.end(); ++eqcit){
              const expression<T> & exp = *eqcit;
              idx2equation_[exp.index].push_back(ep);
            }
        }

      return 0;
    }

    template <typename T>
    int gauss_eliminator<T>::save_equations_mem(
        std::vector<T> & node_vec,
        boost::dynamic_bitset<> & node_flag,
        std::vector<equation<T> > & eq_vec)const
    {
      node_vec = nodes_;
      node_flag = node_flag_;
      eq_vec.clear();
      eq_vec.reserve(es.size());
      for(const_equation_ptr cep = es.begin(); cep != es.end(); ++cep){
          if(cep->state() == 0) continue;
          if(cep->state() == 1){
              // here takes an assumption, this equations is resolved
              continue;
            }
          if(cep->state() == -1){
              std::cerr << "# [error] this equation is conflict " << std::endl;
              std::cerr << *cep << std::endl;
              continue;
            }
          eq_vec.push_back(*cep);
        }
      return 0;
    }

    template <typename T>
    int gauss_eliminator<T>::load_equations_mem(
        const std::vector<T> & node_vec,
        const boost::dynamic_bitset<> & node_flag,
        const std::vector<equation<T> > & eq_vec)
    {
      nodes_ = node_vec;
      node_flag_ = node_flag;

      clear();

      for(size_t eqi = 0; eqi < eq_vec.size(); ++eqi){
          es.push_back(eq_vec[eqi]);

          equation_ptr ep = es.end();
          --ep;
          const equation<T> & eq = *ep;
          prime_idx2equation_[eq.e_vec_.front().index].push_back(ep);

          for(typename equation<T>::eq_const_iterator eqcit = eq.begin();
              eqcit != eq.end(); ++eqcit){
              const expression<T> & exp = *eqcit;
              idx2equation_[exp.index].push_back(ep);
            }
        }
      return 0;
    }

    template <typename T>
    int gauss_eliminator<T>::convert_to_matrix(
        zjucad::matrix::matrix<T> & A,
        zjucad::matrix::matrix<T> & B)const
    {
      std::vector<std::vector<std::pair<size_t, T> > > eqs;
      std::vector<T> b;
      for(typename std::list<equation<T> >::const_iterator lecit = es.begin();
          lecit != es.end(); ++lecit){
          if((*lecit).state() == 0) continue;
          if((*lecit).state() == -1) continue;
          const equation<T> & eq = *lecit;
          std::vector<std::pair<size_t,T> > eq_vec;
          for(typename equation<T>::eq_const_iterator eqcit = eq.begin(); eqcit != eq.end();
              ++eqcit){
              const expression<T> & exp = *eqcit;
              eq_vec.push_back(make_pair(exp.index, exp.coefficient));
            }
          eqs.push_back(eq_vec);
          b.push_back(eq.value());
        }

      A.resize(eqs.size(), nodes_.size());
      B.resize(b.size(),1);

      copy(b.begin(), b.end(), B.begin());
      for(size_t eqi = 0; eqi < eqs.size(); ++eqi){
          const std::vector<std::pair<size_t,T> > & one_eq = eqs[eqi];
          for(size_t exp_i = 0; exp_i < one_eq.size(); ++exp_i){
              A(eqi, one_eq[exp_i].first) = one_eq[exp_i].second;
            }
        }

      return 0;
    }
  }
}
#endif // GAUSS_ELIMINATION_H
