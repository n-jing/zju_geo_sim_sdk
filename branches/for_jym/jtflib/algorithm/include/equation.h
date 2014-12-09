#ifndef JTF_ALGORITHM_EQUATION_H
#define JTF_ALGORITHM_EQUATION_H

#include <map>
#include <vector>
#include <list>
#include <deque>
#include <iostream>
#include <fstream>
#include <boost/unordered_map.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>
#include <zjucad/matrix/matrix.h>

namespace jtf{
  namespace algorithm{


    template <typename T>
    /**
 * @brief store an expression A * i, where A is the coeffieient
 *
 */
    class expression{
    public:
      expression():index(-1),coefficient(0){}
      expression(const size_t & index_, const T & coefficient_)
        : index(index_), coefficient(coefficient_){
        if(static_cast<int>(index_) < 0)
          std::cerr << "# [error] this index is negative." << std::endl;
      }
      size_t index;
      T coefficient;
      /**
   * @brief operator < is used to determin which expression ha
   *
   * @param other input other expression
   * @return bool
   */
      bool operator < (const expression<T> & b) const{
        return index < b.index;
      }

      /**
   * @brief To determin whether coefficent of this expression is zeros
   *
   * @return bool
   */
      bool is_zero() const{ return fabs(coefficient) < 1e-6;}
      int write(std::ostream & os)const;
      int read(std::istream & is);
    };

    template <typename T>
    int expression<T>::write(std::ostream & os)const
    {
      os.write((char *)&index, sizeof(size_t));
      os.write((char *)&coefficient, sizeof(T));
      return 0;
    }

    template <typename T>
    int expression<T>::read(std::istream & is)
    {
      is.read((char *)&index, sizeof(size_t));
      is.read((char *)&coefficient, sizeof(T));
      return 0;
    }

    /**
 * @brief make expression with node_idx and value
 *
 * @param node_idx
 * @param value
 * @return expression<T>
 */
    template <typename T>
    expression<T> make_expression(const size_t & node_idx, const T & value)
    {
      expression<T> temp(node_idx, value);
      return temp;
    }

    /**
 * @brief equation class
 *
 */
    template <typename T>
    class equation{
    public:
      typedef typename std::list<expression<T> >::const_iterator eq_const_iterator;
      typedef typename std::list<expression<T> >::iterator eq_iterator;

      equation():value_(static_cast<T>(0)){}

      eq_const_iterator begin() const {return e_vec_.begin();}
      eq_iterator begin(){return e_vec_.begin();}

      eq_const_iterator end() const {return e_vec_.end();}
      eq_iterator end(){return e_vec_.end();}

      /**
   * @brief used to standardizate the equation, sort the expressions and merge
   *        similar items and normalization to make the first item's coefficient
   *        equal 1
   *
   * @return int
   */
      int standardization(){
        sort_equation();
        merge_similar();
        normalization();
        return 0;
      }

      /**
   * @brief update the equation with given node idx and value
   *
   * @param node_idx input node idx
   * @param node_value input node value
   * @return int if nothing change return 0, or return changes number
   */
      int update(const size_t & node_idx, const  T & node_value);
      /**
   * @brief sort equation accronding to the expressions
   *
   * @return int
   */
      int sort_equation(){
        e_vec_.sort();
        return 0;
      }

      /**
   * @brief merge similar items, WARNING. this function should be used after
   *        sorting.
   *
   * @return int
   */
      int merge_similar();

      /**
   * @brief normalize the equations: to make the first expression's coefficient
   *        equals 1
   *
   * @return int
   */
      int normalization();

      /**
   * @brief get the state of equation:
   *            if there are no expressions,
   *               if value != 0 return -1; error equation
   *               else return 0; cleared
   *            else
   *               if there is one expression return 1; calculated finial node
   *               else return 2; not finished
   *
   * @return int
   */
      int state() const;


      /**
   * @brief define the minus operation
   *
   * @param b
   * @return equation<T>
   */
      equation<T> & operator -= (const equation<T> & b);

      equation<T> & operator *= (const T& b);

      /**
     * @brief define the output operation
     *
     * @param eq
     * @return equation<T>
     */
      friend std::ostream& operator << (std::ostream &output,
                                        const equation<T> &eq)
      {
        if(eq.e_vec_.empty()){
            output << "# ------ empty expressions with value = " << eq.value() << std::endl;
          }else{
            output << "# ------ SUM coeff * index , val" << std::endl
                   << "# ------ ";
            for(equation<T>::eq_const_iterator eqcit = eq.begin(); eqcit != eq.end();){
                const expression<T> & exp = *eqcit;
                output << exp.coefficient << "*X" << exp.index;
                ++eqcit;
                if(eqcit == eq.end())
                  output << " = ";
                else
                  output << " + ";
              }
            output << eq.value() << std::endl;
          }
        return output;
      }

      int write(std::ostream & ofs)const;
      int read(std::istream & ifs);

      /**
   * @brief get the first expression idx
   *
   * @return size_t
   */
      size_t get_prime_idx() const {
        if(e_vec_.empty())
          return -1;
        else
          return e_vec_.front().index;
      }

      int get_coefficient(const size_t & index, T& coeff)const {
        for(typename std::list<expression<T> >::const_iterator lcit = e_vec_.begin();
            lcit != e_vec_.end(); ++lcit){
            if(lcit->index == index){
                coeff = lcit->coefficient;
                return 0;
              }
          }
        return __LINE__;
      }

      int add_expression(const expression<T> & exp){
        e_vec_.push_back(exp);
        return 0;
      }

      T& value() {return value_;}
      const T& value() const {return value_;}

      std::list<expression<T> > e_vec_;
    private:
      T value_;
    };

    template <typename T>
    int equation<T>::write(std::ostream & ofs)const
    {
      const size_t item_num = e_vec_.size();
      ofs.write((char*)&item_num, sizeof(size_t));
      ofs.write((char*)&value_, sizeof(T));
      for(eq_const_iterator expit = e_vec_.begin();
          expit != e_vec_.end(); ++expit){
          expit->write(ofs);
        }
      return 0;
    }

    template <typename T>
    int equation<T>::read(std::istream & ifs)
    {
      size_t item_num = 0;
      ifs.read((char*)&item_num, sizeof(size_t));
      ifs.read((char*)&value_, sizeof(T));
      expression<T> exp;
      for(size_t expi = 0; expi < item_num; ++expi){
          exp.read(ifs);
          add_expression(exp);
        }
      return 0;
    }


    template <typename T>
    int equation<T>::merge_similar(){
      typedef typename std::list<expression<T> >::iterator leit;
      leit current = e_vec_.begin();
      leit next = current;

      //  while(current != e_vec_.end()){
      //    next = current;
      //    ++next;
      //    if(next == e_vec_.end()) break;
      //    if(next->index > current->index){
      //      current = next;
      //      continue;
      //    }else if(next->index == current->index){
      //      current->coefficient += next->coefficient;
      //      e_vec_.erase(next);
      //      if(fabs(current->coefficient) < 1e-6){
      //        e_vec_.erase(current++);
      //      }
      //    }else{
      //      std::cerr << "# [error] merge similar function should only be called "
      //                << "after sorting." << std::endl;
      //      return __LINE__;
      //    }
      //  }
      ++next;
      while(next != e_vec_.end()){
          if(fabs(current->coefficient) < 1e-6){
              e_vec_.erase(current++);
              next=current;
              ++next;
              continue;
            }
          if(next->index > current->index) {
              current = next++;
              continue;
            }else if(next->index == current->index){
              current->coefficient += next->coefficient;
              e_vec_.erase(next++);
              if(fabs(current->coefficient) < 1e-6){
                  e_vec_.erase(current++);
                  if(current == e_vec_.end()) break;
                  ++next;
                }
            }else{
              std::cerr << "# [error] merge similar function should only be called "
                        << "after sorting." << std::endl;
              return __LINE__;
            }
        }
      return 0;
    }

    template <typename T>
    int equation<T>::normalization()
    {
      const T coeff = e_vec_.front().coefficient;
      if(fabs(coeff) < 1e-6){
          if(e_vec_.empty()){
              //std::cerr << "# [info] this equation is empty." << std::endl;
              return 0;
            }else
            std::cerr << "# [error] this expression should be removed." << std::endl;
          return __LINE__;
        }

      value() /= coeff;
      for(typename std::list<expression<T> >::iterator lit = e_vec_.begin();
          lit != e_vec_.end(); ++lit){
          expression<T> & ep = *lit;
          ep.coefficient /= coeff;
        }
      return 0;
    }

    template <typename T>
    int equation<T>::state() const
    {
      if(e_vec_.empty()){
          if(fabs(value()) < 1e-8)
            return 0; // is cleared
          return -1; // is conflicted
        }else{
          if(e_vec_.size() == 1){
              if(fabs(e_vec_.front().coefficient)<1e-8 && fabs(value())<1e-8)
                return 0;
              if(fabs(e_vec_.front().coefficient)<1e-8 && fabs(value())>1e-8)
                return -1;
              return 1; // finial variant
            }
          else
            return 2; // not ready
        }
    }

    template <typename T>
    equation<T> & equation<T>::operator *= (const T& b)
    {
      for(typename std::list<expression<T> >::iterator leit = e_vec_.begin();
          leit != e_vec_.end(); ++leit){
          leit->coefficient *= b;
        }
      value() *= b;
      return *this;
    }

    template <typename T>
    equation<T> & equation<T>::operator -= (const equation<T> & b )
    {
      if(&b == this) {
          e_vec_.clear();
          value() = 0;
          return *this;
        }

      //  for(typename std::list<expression<T> >::const_iterator lecit_b =
      //      b.e_vec_.begin(); lecit_b != b.e_vec_.end(); ++lecit_b){
      //    expression<T> exp = *lecit_b;
      //    exp.coefficient *= -1;
      //    add_expression(exp);
      //  }
      for(typename std::list<expression<T> >::const_iterator lecit_b =
          b.e_vec_.begin(); lecit_b != b.e_vec_.end(); ++lecit_b){
          const expression<T> & exp = *lecit_b;
          const size_t &node_idx = exp.index;
           if(fabs(exp.coefficient) < 1e-6) continue;
          bool is_found = false;

          for(typename std::list<expression<T> >::iterator leit_a = e_vec_.begin();
              leit_a != e_vec_.end(); ++leit_a){
              expression<T> & exp_a = *leit_a;
              if(exp_a.index == node_idx){
                  exp_a.coefficient -= exp.coefficient;
                  is_found = true;
                  if(fabs(exp_a.coefficient) < 1e-6) {
                      e_vec_.erase(leit_a);
                      break;
                    }
                }
            }
          if(!is_found){
              e_vec_.push_back(make_expression(node_idx, -1 * exp.coefficient));
            }
        }

      value() -= b.value();
      standardization();
      //  sort_equation();
      //  normalization();

      return *this;
    }

    template <typename T>
    int equation<T>::update(const size_t & node_idx, const  T & node_value)
    {
      int changes = 0;
      for(typename std::list<expression<T> >::iterator leit = e_vec_.begin();
          leit != e_vec_.end();){
          expression<T> & exp = *leit;
          if(exp.index == node_idx){
              value() -= exp.coefficient * node_value;
              e_vec_.erase(leit++);
              ++changes;
            }
          ++leit;
        }

      standardization();
      return changes;
    }
  }
}

#endif
