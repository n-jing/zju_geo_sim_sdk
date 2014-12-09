#ifndef EDGE_SXX_H
#define EDGE_SXX_H

#include <zjucad/matrix/matrix.h>
#include <cassert>
#include <vector>
#include <set>

#include "basic_define.h"
#include "Property.h"
//class of edge
namespace sxx
{
  class Edge
  {
  public:
    Edge():edge_(std::make_pair(0, 0)), mid_index_(0), size_field_(0.0)
    {
    }

    Edge(const std::pair<size_t, size_t> &e):edge_(e), mid_index_(0), size_field_(0.0)
    {
      assert(e.first < e.second);
    }

    Edge(const size_t id1, const size_t id2)
    {
      assert(id1 != id2);
      if(id1 < id2)
        {
          edge_.first = id1;
          edge_.second = id2;
        }
      else
        {
          edge_.first = id2;
          edge_.second = id1;
        }
      mid_index_ = 0;
      size_field_ = 0.0;
    }

    Edge(const Edge &e)
    {
      edge_ = e.get_edge();
      mid_index_ = e.get_mid_index();
      size_field_ = e.get_size_field();
    }

    const std::pair<size_t, size_t> &get_edge() const
    {
      return edge_;
    }

    std::pair<size_t, size_t> &get_edge()
    {
      return edge_;
    }

    bool operator==(const Edge &e)
    {
      return edge_ == e.get_edge();
    }

    void set_length(const double length)
    {
      length_ = length;
    }

    const double &get_length() const
    {
      return length_;
    }

    const edgeframe_list_type &get_edgeframe_list() const
    {
      return edgeframe_list_;
    }

    edgeframe_list_type &get_edgeframe_list()
    {
      return edgeframe_list_;
    }

    const UV_property &get_uv_property() const
    {
      return uv_property_;
    }

    UV_property &get_uv_property()
    {
      return uv_property_;
    }

    size_t get_mid_index() const
    {
      return mid_index_;
    }

    void set_mid_index(size_t mid)
    {
      mid_index_ = mid;
    }


    double get_size_field()
    {
      return size_field_;
    }

    void set_size_field(double size_field)
    {
      size_field_ = size_field;
    }

    void set_size_field()
    {
      size_field_ = length_ * (zjucad::matrix::norm(uv_property_.get_uv()));
    }

    virtual ~Edge()
    {
    }

  private:
    std::pair<size_t, size_t> edge_;  //the first is less than the second
    edgeframe_list_type edgeframe_list_;
    UV_property uv_property_;
    double length_;
    size_t mid_index_;
    double size_field_;
  };
}

#endif // EDGE_SXX_H
