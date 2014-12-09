#ifndef BASIC_DEFINE_SXX_H
#define BASIC_DEFINE_SXX_H

#include <iostream>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <zjucad/matrix/matrix.h>
#include <cmath>
#include <iterator>
#include <list>
#include <set>

//#include "Vertex.h"
//#include "Edge.h"
//#include "Tri_face.h"
//#include "Property.h"

namespace sxx
{
  class Vertex;
  class Edge;
  class Tri_face;
  class Property;
  class Frame_property;
//  size_t hash_value(const Edge &edge);

//  size_t hash_value(const Tri_face &face);

  int load_obj(const char * file_name,
                    std::vector<std::vector<double> > &points,
                    std::vector<size_t> &faces);

  template<typename T>
  bool is_element_in_vec(const T &element, const std::vector<T> &vec)
  {
    for(size_t i = 0; i < vec.size(); ++i)
      if(element == vec[i])
        return true;
    return false;
  }

  typedef std::list<Edge> edgelist_type;
  typedef std::list<Edge>::iterator edgeit_type;
  typedef std::list<Tri_face> trifacelist_type;
  typedef std::list<Tri_face>::iterator trifaceit_type;
  typedef std::list<Frame_property> framelist_type;
  typedef std::list<Frame_property>::iterator frameit_type;
  typedef boost::unordered_map<std::pair<size_t, size_t>, edgeit_type> edge2it_type;
  typedef std::pair<trifaceit_type, Frame_property> edgeframe_type;
  typedef std::list<edgeframe_type> edgeframe_list_type;
}
#endif // BASIC_DEFINE_SXX_H
