#ifndef TRI_FACE_SXX_H
#define TRI_FACE_SXX_H

#include <zjucad/matrix/matrix.h>
#include <cassert>
#include <vector>
#include <set>

#include "Property.h"
#include "basic_define.h"

//class of triangle_face
namespace sxx
{
  class Tri_face
  {
  public:
    Tri_face(const size_t id1, const size_t id2, const size_t id3):face_(3, 0), face_index_(-1)
    {
      assert(id1 != id2 && id2 != id3 && id1 != id3);
      face_[0] = id1;
      face_[1] = id2;
      face_[2] = id3;
    }

    Tri_face(const std::vector<size_t> &face):face_index_(-1)
    {
      assert(face.size() == 3);
      face_ = face;
    }

    Tri_face(const Tri_face &tri_face)
    {
      face_ = tri_face.get_face();
      face_index_ = tri_face.get_face_index();
    }

    Tri_face(const zjucad::matrix::matrix<size_t> &face):face_(3, 0), face_index_(-1)
    {
      assert(face.size() == 3);
      for(size_t i = 0; i < face.size(); ++i)
        face_[i] = face[i];
    }

    const std::vector<size_t> &get_face() const
    {
      return face_;
    }

    std::vector<size_t> &get_face()
    {
      return face_;
    }

    const std::vector<edgeit_type> &get_adjacent_edges() const
    {
      return adjacent_edges_;
    }

    std::vector<edgeit_type> &get_adjacent_edges()
    {
      return adjacent_edges_;
    }

    void set_face_index(const size_t num)
    {
      face_index_ = num;
    }

    const int get_face_index() const
    {
      return face_index_;
    }

    bool operator==(const Tri_face &f)
    {
      return face_ == f.get_face();
    }

    virtual ~Tri_face()
    {
    }

  private:
    std::vector<size_t> face_;
    std::vector<edgeit_type> adjacent_edges_;
    int face_index_;  // it is used when write info to evfile
  };
}
#endif // TRI_FACE_SXX_H
