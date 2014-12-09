#include "../include/gauss_elimination.h"
#include <iostream>

using namespace std;
using namespace boost;

int main()
{
  vector<double> x(5);
  boost::dynamic_bitset<> flag(5);
  jtf::algorithm::gauss_eliminator<double> ge(x,flag);

  {jtf::algorithm::equation<double> eqn;
    eqn.add_expression(jtf::algorithm::make_expression(1,1.0));
    eqn.add_expression(jtf::algorithm::make_expression(3,1.0));
    ge.add_equation(eqn);
  }
  {
    jtf::algorithm::equation<double> eqn;
    eqn.add_expression(jtf::algorithm::make_expression(3,-1.0));
    eqn.add_expression(jtf::algorithm::make_expression(4,-1.0));
    ge.add_equation(eqn);
  }

  {
    jtf::algorithm::equation<double> eqn;
    eqn.add_expression(jtf::algorithm::make_expression(2,-2.0));
    eqn.add_expression(jtf::algorithm::make_expression(3,-1.0));
    ge.add_equation(eqn);
  }

  {
    jtf::algorithm::equation<double> eqn;
    eqn.add_expression(jtf::algorithm::make_expression(4,-2.0));
    ge.add_equation(eqn);
  }


  for(const auto & one_eqn : ge){
      cerr << "eqn" << endl;
      for(const auto & one_exp : one_eqn){
          cerr << one_exp.index << " ";
        }
      cerr << endl;
      for(const auto & one_exp : one_eqn){
          cerr << one_exp.coefficient << " ";
        }
      cerr << endl;
    }

  for(size_t i = 0; i < flag.size(); ++i){
      if(flag[i])
        cerr << i << " = " << x[i] << endl;
    }
  return 0;
}
