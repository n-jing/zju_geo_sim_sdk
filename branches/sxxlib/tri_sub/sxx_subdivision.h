#ifndef MYSUBDIVISOIN_SXX_H
#define MYSUBDIVISOIN_SXX_H

#include <cstdio>
#include <fstream>
#include "basic_define.h"
#include "Vertex.h"
#include "Edge.h"
#include "Tri_face.h"

namespace sxx
{
  size_t hash_value(const Edge &edge);

  size_t hash_value(const Tri_face &face);

  const double COS45_  = 0.70710678118655;

  class Subdivision
  {
  public:
    Subdivision():threshhold_(0.0)
    {
      vertex_vec_.clear();
      edge_list_.clear();
      triface_list_.clear();
      edge2it_.clear();
      feature_line_.clear();
      evfile_flag_ = false;
      esfile_flag_ = false;
      flfile_flag_ = false;
    }


    int import_mesh(const char *obj_file, const char *es_file = NULL,
                    const char *ev_file = NULL, const char *fl_file = NULL);

    int import_mesh(const zjucad::matrix::matrix<double> &points,
                    const zjucad::matrix::matrix<size_t> &tri_face);

    int export_mesh(zjucad::matrix::matrix<double> &points,
                    zjucad::matrix::matrix<size_t> &tri_face);

    int export_mesh(const char *obj_file, const char *es_file = NULL,
                    const char *ev_file = NULL, const char *fl_file = NULL);

    int subdived_edges(const std::vector<std::pair<size_t, size_t> > &sub_edges, std::vector<size_t> &insert_points);

    const double get_threshhold() const
    {
      return threshhold_;
    }

    void set_threshhold(const double threshhold)
    {
      threshhold_ = threshhold;
    }

    int adaptive_subdivided();

    int flip_edges(const std::vector<std::pair<size_t, size_t> > &edges);

    virtual ~Subdivision()
    {
      vertex_vec_.clear();
      edge_list_.clear();
      triface_list_.clear();
      edge2it_.clear();
      feature_line_.clear();
    }

  private:
    int load_esfile(const char *es_file);

    int save_esfile(const char *es_file);

    int load_evfile(const char *ev_file);

    int save_evfile(const char *ev_file);

    int load_flfile(const char *fl_file);

    int save_flfile(const char *fl_file);

    int save_objfile(const char *obj_file);

    int add_face(const size_t id0, const size_t id1, const size_t id2, trifaceit_type &cur_face_it);

    int add_edge(const size_t id0, const size_t id1, const trifaceit_type &tri_face_it, edgeit_type &cur_edge_it);

    int delete_face(const trifaceit_type &triface_it);

    int del_face_from_edge(const edgeit_type &edge_it, const trifaceit_type &triface_it);

    int insert_uv_property(const size_t id0, const size_t id1,
                           const zjucad::matrix::matrix<double> &uv_value);

    int insert_frame_property(const size_t id0, const size_t id1, const size_t face_no,
                              const zjucad::matrix::matrix<double> &value);

    int align_frame(const zjucad::matrix::matrix<double> &base_vector,
                    const zjucad::matrix::matrix<double> &normal_vector,
                    zjucad::matrix::matrix<double> &aligned_frame);

    int get_face_three_frames(const trifaceit_type &triface_it,
                              zjucad::matrix::matrix<double> &face_frames);

    int get_face_uv_values(const trifaceit_type &triface_it,
                           zjucad::matrix::matrix<double> &uv_values);

    int get_edges_subdivided(const trifaceit_type &triface_it,
                             std::vector<bool> &edge_flag);

    int sub_one_edge(const trifaceit_type &triface_it, const std::vector<bool> &edge_flag);

    int sub_two_edge(const trifaceit_type &triface_it, const std::vector<bool> &edge_flag);

    int sub_three_edge(const trifaceit_type &triface_it, const std::vector<bool> &edge_flag);

    int adjust_points_index(const trifaceit_type &triface_it, const std::vector<bool> &edge_flag);

    int update_adjacent_edges(const std::vector<size_t> &face,
                              std::vector<edgeit_type> &adjacent_edges);

    int insert_mid_point(const edgeit_type &edge_it);

    int insert_new_edge(const size_t id0, const size_t id1,
                        const bool flag, const trifaceit_type &cur_face_it,
                        const trifaceit_type &old_face_it,
                        const zjucad::matrix::matrix<double> &uv_value, const zjucad::matrix::matrix<double> &frame_value,
                        const double *barycentric, edgeit_type &cur_edge_it);

    int insert_new_face(const trifaceit_type &triface_it, const size_t id0, const size_t id1,
                        const size_t id2, const size_t sub_edge_num, const double *barycentric,
                        const zjucad::matrix::matrix<double> &uv_value, const zjucad::matrix::matrix<double> &frame_value);

    int update_feature_line();

    int delete_remain_edges();

    int edge_flip(const edgeit_type &edge_it);

    std::vector<Vertex> vertex_vec_;
    edgelist_type edge_list_;
    trifacelist_type triface_list_;
    edge2it_type edge2it_;
    bool esfile_flag_, evfile_flag_, flfile_flag_;
    double threshhold_;
    std::vector<std::list<size_t> > feature_line_;
  };
}
#endif // MYSUBDIVISOIN_SXX_H
