#ifndef VERTEX_SXX_H
#define VERTEX_SXX_H

#include <zjucad/matrix/matrix.h>
#include <iostream>
#include <vector>
#include <cassert>

//class of vertex
namespace sxx
{
  class Vertex
  {
  public:
    Vertex():coord_(3, 1)
    {
      coord_[0] = 0.0;
      coord_[1] = 0.0;
      coord_[2] = 0.0;
    }

    Vertex(const double x, const double y, const double z):coord_(3, 1)
    {
      coord_[0] = x;
      coord_[1] = y;
      coord_[2] = z;
    }

    Vertex(const std::vector<double> &ver):coord_(3, 1)
    {
      coord_[0] = ver[0];
      coord_[1] = ver[1];
      coord_[2] = ver[2];
    }

    Vertex(const zjucad::matrix::matrix<double> &coord):coord_(coord)
    {
      assert(coord.size() == 3);
    }

    Vertex(const Vertex& ver)
    {
      this->coord_ = ver.get_coord();
    }

    inline zjucad::matrix::matrix<double> &get_coord()
    {
      return coord_;
    }

    inline const zjucad::matrix::matrix<double> &get_coord() const
    {
      return coord_;
    }

    void inline set_coord(const zjucad::matrix::matrix<double> &coord)
    {
      assert(coord.size() == 3);
      coord_ = coord;
    }

    virtual ~Vertex()
    {
    }

    private:
      zjucad::matrix::matrix<double> coord_;

  };
}

#endif // VERTEX_SXX_H
