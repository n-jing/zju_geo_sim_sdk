#include <string>
#include <stdio.h>
#include <set>
#include <iostream>
#include <fstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/shared_ptr.hpp>

#include "../include/io.h"
#include <map>

#include "mat_data_tree.h"

using namespace std;
using boost::property_tree::ptree;

using namespace zjucad::matrix;
using namespace hj::data_tree;
#define FL , __FILE__, __LINE__

int main(int argc, char *argv[])
{
  const string hello = "hello", world = "world";
  data_tree dt;
  dt.put("str-0", hello);
  dt.put("str-1", world);
  float rgb[] = {50, 24, 190};
  dt.put("rgb/r", rgb[0]);
  dt.put("rgb/g", rgb[1]);
  dt.put("rgb/b", rgb[2]);

  type_op_map prt_tr = default_print_trait();
  print_writer pw(prt_tr);
  pw.write(cout, dt);

  cout << *dt.get<float>("rgb/b") << endl;
  return 0;
}
