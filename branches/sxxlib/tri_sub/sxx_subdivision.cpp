#include "sxx_subdivision.h"

using namespace sxx;
using namespace std;

size_t sxx::hash_value(const Edge &edge)
{
  boost::hash<std::pair<size_t, size_t> > t_hash;
  return t_hash(edge.get_edge());
}

size_t sxx::hash_value(const Tri_face &face)
{
  boost::hash<std::vector<size_t> > t_hash;
  return t_hash(face.get_face());
}

int Subdivision::import_mesh(const char *obj_file, const char *es_file = NULL,
                             const char *ev_file = NULL, const char *fl_file = NULL)
{
  vector<vector<double> > points;
  vector<size_t> faces;
  if(load_obj(obj_file, points, faces))
    {
      cerr << "can not open obj file" << endl;
      esfile_flag_ = false;
      evfile_flag_ = false;
      flfile_flag_ = false;
      return 1;
    }
  assert(faces.size() % 3 == 0);
  for(size_t i = 0;i < points.size(); ++i)
    vertex_vec_.push_back(Vertex(points[i]));
  trifaceit_type cur_face_it;
  for(size_t i = 0; i < faces.size(); i += 3)
    {
      add_face(faces[i], faces[i+1], faces[i+2], cur_face_it);
      cur_face_it->set_face_index(i / 3);
    }
  if(!load_esfile(es_file))
    esfile_flag_ = true;
  if(!load_evfile(ev_file))
    evfile_flag_ = true;
  if(!load_flfile(fl_file))
    flfile_flag_ = true;
  if(esfile_flag_)
    {
      for(edgeit_type edge_it = edge_list_.begin();
          edge_it != edge_list_.end(); ++edge_it)
        edge_it->set_size_field();
    }
  return 0;
}

int Subdivision::adaptive_subdivided()
{
  vector<bool> edge_flag(3, false);
  int sub_edge_num;
  trifaceit_type triface_it, triface_it1;
  size_t final_face_index = 0;
  for(triface_it = triface_list_.begin();
      triface_it != triface_list_.end();)
    {
      sub_edge_num = get_edges_subdivided(triface_it, edge_flag);
      if(sub_edge_num == 0)
        {
          triface_it->set_face_index(final_face_index);
          ++final_face_index;
          ++triface_it;
        }
      else if(sub_edge_num == 1)
        {
          sub_one_edge(triface_it, edge_flag);
          triface_it1 = triface_it;
          ++triface_it;
          delete_face(triface_it1);
        }
      else if(sub_edge_num == 2)
        {
          sub_two_edge(triface_it, edge_flag);
          triface_it1 = triface_it;
          ++triface_it;
          delete_face(triface_it1);
        }
      else if(sub_edge_num == 3)
        {
          sub_three_edge(triface_it, edge_flag);
          triface_it1 = triface_it;
          ++triface_it;
          delete_face(triface_it1);
        }
      else
        cerr << "error" << endl;
    }
//  if(flfile_flag_)
//    {
//      update_feature_line();
//      delete_remain_edges();
//    }
  return 0;
}

int Subdivision::export_mesh(const char *obj_file, const char *es_file = NULL,
                             const char *ev_file = NULL, const char *fl_file = NULL)
{
  save_objfile(obj_file);
  if(esfile_flag_)
    save_esfile(es_file);
  if(evfile_flag_)
    save_evfile(ev_file);
  if(flfile_flag_)
    save_flfile(fl_file);
  return 0;
}

int Subdivision::add_face(const size_t id0, const size_t id1, const size_t id2, trifaceit_type &cur_face_it)
{
  triface_list_.push_back(Tri_face(id0, id1, id2));
  cur_face_it = --triface_list_.end();
  vector<edgeit_type> &adjacent_edges = cur_face_it->get_adjacent_edges();
  edgeit_type cur_edge_it;
  vector<size_t> temp_vec(3, 0);
  temp_vec[0] = id0;
  temp_vec[1] = id1;
  temp_vec[2] = id2;
  for(size_t i = 0; i < 3; ++i)
    {
      add_edge(temp_vec[i], temp_vec[(i+1)%3], cur_face_it, cur_edge_it);
      adjacent_edges.push_back(cur_edge_it);
    }
  assert(adjacent_edges.size() == 3);
//  cout << "add face " << id0 << " " << id1 << " " << id2 << " OK" << endl;
  return 0;
}

int Subdivision::add_edge(const size_t id0, const size_t id1, const trifaceit_type &tri_face_it, edgeit_type &cur_edge_it)
{
  assert(id0 != id1);
  pair<size_t, size_t> edge_pair(id0, id1);
  if(id0 > id1)
    swap(edge_pair.first, edge_pair.second);
  const edge2it_type::iterator &edge2it_it = edge2it_.find(edge_pair);
  if(edge2it_it == edge2it_.end())
    {
      edge_list_.push_back(Edge(edge_pair));
      cur_edge_it = --edge_list_.end();
      edge2it_.insert(make_pair(edge_pair, cur_edge_it));
      cur_edge_it->set_length(zjucad::matrix::norm(vertex_vec_[id0].get_coord() -
                                                   vertex_vec_[id1].get_coord()));
      edgeframe_list_type &edgeframe_list = cur_edge_it->get_edgeframe_list();
      edgeframe_list.push_back(make_pair(tri_face_it, Frame_property()));
    }
  else
    {
      cur_edge_it = edge2it_it->second;
      edgeframe_list_type &edgeframe_list = cur_edge_it->get_edgeframe_list();
      edgeframe_list.push_back(make_pair(tri_face_it, Frame_property()));
      assert(edgeframe_list.size() == 2);
    }
  return 0;
}

int Subdivision::delete_face(const trifaceit_type &triface_it)
{
  const vector<edgeit_type> &adjacent_edges = triface_it->get_adjacent_edges();
  for(size_t i = 0; i < adjacent_edges.size(); ++i)
    {
      del_face_from_edge(adjacent_edges[i], triface_it);
    }
  triface_list_.erase(triface_it);
  return 0;
}

int Subdivision::del_face_from_edge(const edgeit_type &edge_it, const trifaceit_type &triface_it)
{
  edgeframe_list_type &edgeframe_list = edge_it->get_edgeframe_list();
  for(edgeframe_list_type::iterator edgeframe_it = edgeframe_list.begin();
      edgeframe_it != edgeframe_list.end(); ++edgeframe_it)
    {
      if(edgeframe_it->first->get_face() == triface_it->get_face())
        {
          edgeframe_list.erase(edgeframe_it);
          if(edgeframe_list.size() == 0 && !flfile_flag_)
            {
              edge_list_.erase(edge_it);
              edge2it_.erase(edge_it->get_edge());
            }
          break;
        }
    }
  return 0;
}


int Subdivision::load_esfile(const char *es_file)
{
  if(!es_file)
    return 2;
  ifstream ifs(es_file);
  if( !ifs )
    {
      cerr << "load es file error!" << endl;
      return 1;
    }
  size_t num, id0, id1, trash;
  zjucad::matrix::matrix<double> uv(2, 1);
  ifs >> num;
  for(size_t i = 0; i < num; ++i)
    {
      ifs >> trash >> id0 >> id1 >> uv[0] >> uv[1];
      if(id0 > id1)
        swap(id0, id1);
      insert_uv_property(id0, id1, uv);
    }
  ifs.close();
  cout << "load es file end!" << endl;
  return 0;
}

int Subdivision::insert_uv_property(const size_t id0, const size_t id1,
                                    const zjucad::matrix::matrix<double> &uv_value)
{
  assert(id0 < id1);
  edge2it_type::iterator edge2it_it = edge2it_.find(make_pair(id0, id1));
  assert(edge2it_it != edge2it_.end());
  (edge2it_it->second->get_uv_property()).set_uv_value(uv_value);
}

int Subdivision::load_evfile(const char *ev_file)
{
  if(!ev_file)
    return 2;
  ifstream ifs(ev_file);
  if(!ifs)
    {
      cerr << "load ev file error!" << endl;
      return 1;
    }
  size_t num, trash, id0, id1, face_num, face_no;
  zjucad::matrix::matrix<double> frame_value(3, 1);
  ifs >> num;
  for(size_t i = 0; i < num; ++i)
    {
      ifs >> trash >> id0 >> id1 >> face_num;
      if(id0 > id1)
        swap(id0, id1);
      for(size_t j = 0; j < face_num; ++j)
        {
          ifs >> face_no >> frame_value[0] >> frame_value[1] >> frame_value[2];
          insert_frame_property(id0, id1, face_no, frame_value);
        }
    }
  ifs.close();
  cout << "load ev file end!" << endl;
  return 0;
}

int Subdivision::insert_frame_property(const size_t id0, const size_t id1, const size_t face_no,
                                       const zjucad::matrix::matrix<double> &value)
{
  assert(id0 < id1);
  edge2it_type::iterator edge2it_it = edge2it_.find(make_pair(id0, id1));
  assert(edge2it_it != edge2it_.end());
  edgeframe_list_type &edgeframe_list = edge2it_it->second->get_edgeframe_list();
  for(edgeframe_list_type::iterator edgeframe_it = edgeframe_list.begin();
      edgeframe_it != edgeframe_list.end(); ++edgeframe_it)
    {
      if((edgeframe_it->first)->get_face_index() == face_no)
        {
          (edgeframe_it->second).set_frame(value);
          break;
        }
    }
  return 0;
}

int Subdivision::save_objfile(const char *obj_file)
{
  ofstream ofs(obj_file);
  if(!ofs)
    {
      cerr << "save obj file error!";
      return 1;
    }
  for(size_t i = 0; i < vertex_vec_.size(); ++i)
    {
      const zjucad::matrix::matrix<double> &ver = vertex_vec_[i].get_coord();
      ofs << "v " << ver[0] << " "
                  << ver[1] << " "
                  << ver[2] << endl;
    }

  for(trifaceit_type triface_it = triface_list_.begin();
      triface_it != triface_list_.end(); ++triface_it)
    {
      const vector<size_t> &face = triface_it->get_face();
      ofs << "f " << face[0] + 1 << " "
                  << face[1] + 1 << " "
                  << face[2] + 1 << endl;
    }
  cout << "save obj file end" << endl;
  return 0;
}

int Subdivision::save_esfile(const char *es_file)
{
  ofstream ofs(es_file);
  if(!ofs)
    {
      cerr << "save es file error!";
      return 1;
    }
  size_t cnt = 0;
  ofs << edge_list_.size() << endl;
  for(edgeit_type edge_it = edge_list_.begin();
      edge_it != edge_list_.end(); ++edge_it)
    {
      const pair<size_t, size_t> &edge_pair = edge_it->get_edge();
      const zjucad::matrix::matrix<double> &uv_value = (edge_it->get_uv_property()).get_uv();
      ofs << cnt << " " << edge_pair.first << " " << edge_pair.second
          << " " << uv_value[0] << " " << uv_value[1] << endl;
      ++cnt;

    }
  cout << "save es file end!" << endl;
  return 0;
}

int Subdivision::save_evfile(const char *ev_file)
{
  ofstream ofs(ev_file);
  if(!ofs)
    {
      cerr << "save ev file error!";
      return 1;
    }
  size_t cnt = 0;
  ofs << edge_list_.size() << endl;
  for(edgeit_type edge_it = edge_list_.begin();
      edge_it != edge_list_.end(); ++edge_it)
    {
      const pair<size_t, size_t> &edge_pair = edge_it->get_edge();
      const edgeframe_list_type &edgeframe_list = edge_it->get_edgeframe_list();
      assert(edgeframe_list.size() > 0 && edgeframe_list.size() < 3);
      ofs << cnt << " " << edge_pair.first << " " << edge_pair.second << " " << edgeframe_list.size() << endl;
      for(edgeframe_list_type::const_iterator edgeframe_it = edgeframe_list.begin();
          edgeframe_it != edgeframe_list.end(); ++edgeframe_it)
        {
          const zjucad::matrix::matrix<double> &frame = (edgeframe_it->second).get_frame();
          ofs << edgeframe_it->first->get_face_index() << " "  << frame[0]
              << " " << frame[1] << " " << frame[2] << endl;
        }
      ++cnt;
    }
  cout << "save ev file end!" << endl;
  return 0;
}

int Subdivision::align_frame(const zjucad::matrix::matrix<double> &base_vector, // identity vector
                             const zjucad::matrix::matrix<double> &normal_vector, // identity vector
                             zjucad::matrix::matrix<double> &aligned_frame)

{
  assert(base_vector.size() == 3 && normal_vector.size() == 3 && aligned_frame.size() == 3);
  double len = zjucad::matrix::norm(aligned_frame);
  assert(len >= 1e-6);
  aligned_frame /= len;
  double tempcos = zjucad::matrix::dot(base_vector, aligned_frame);
  if( abs(tempcos)>= COS45_ )
    {
      if(tempcos < 0)
        aligned_frame = -aligned_frame;
    }
  else
    {
      zjucad::matrix::matrix<double> roate_frame = zjucad::matrix::cross(normal_vector, aligned_frame);
      double len = zjucad::matrix::norm(roate_frame);
      assert(len >= 1e-6);
      roate_frame /= len;
      double temp_cross = zjucad::matrix::dot(base_vector, roate_frame);
      if(temp_cross >= 0)
        {
          aligned_frame = roate_frame;
        }
      else
        {
          aligned_frame = -roate_frame;
        }
    }
  return 0;
}

int Subdivision::get_face_three_frames(const trifaceit_type &triface_it,
                                       zjucad::matrix::matrix<double> &face_frames)
{
  assert(face_frames.size(1) == 3 && face_frames.size(2) ==3);
  zjucad::matrix::matrix<double> frame(3, 1);
  const vector<edgeit_type> &adjacent_edges = triface_it->get_adjacent_edges();
  const vector<size_t> &face = triface_it->get_face();
  zjucad::matrix::matrix<double> base_vector = vertex_vec_[face[1]].get_coord() - vertex_vec_[face[0]].get_coord();
  zjucad::matrix::matrix<double> temp_vector = vertex_vec_[face[2]].get_coord() - vertex_vec_[face[0]].get_coord();

  double len = zjucad::matrix::norm(base_vector);
  assert(len >= 1e-6);
  base_vector /= len;

  len = zjucad::matrix::norm(temp_vector);
  assert(len >= 1e-6);
  temp_vector /= len;

  zjucad::matrix::matrix<double> normal_vector = zjucad::matrix::cross(base_vector, temp_vector);
  len = zjucad::matrix::norm(normal_vector);
    assert(len >= 1e-6);
  normal_vector /= len;

  assert(adjacent_edges.size() == 3);
  size_t cnt = 0;
  vector<pair<size_t, size_t> > edge_vec(3);
  for(size_t i = 0; i < 3 ; ++i)
    edge_vec[i] = adjacent_edges[i]->get_edge();
  for(vector<edgeit_type>::const_iterator adj_edge_it = adjacent_edges.begin();
      adj_edge_it != adjacent_edges.end(); ++adj_edge_it, ++cnt)
    {
      const edgeframe_list_type &edgeframe_list = (*adj_edge_it)->get_edgeframe_list();
      size_t num = 0;
      for(edgeframe_list_type::const_iterator edgeframe_it = edgeframe_list.begin();
          edgeframe_it != edgeframe_list.end(); ++edgeframe_it, ++num)
        {
          if(edgeframe_it->first->get_face() == triface_it->get_face())
            {
              frame = (edgeframe_it->second).get_frame();
              break;
            }
        }
      align_frame(base_vector, normal_vector, frame);
      face_frames(zjucad::matrix::colon(), cnt) = frame;
    }
  return 0;
}

int Subdivision::get_face_uv_values(const trifaceit_type &triface_it,
                                    zjucad::matrix::matrix<double> &uv_values)
{
  assert(uv_values.size(1) == 2 && uv_values.size(2) == 3);
  const vector<edgeit_type> &adjacent_edges = triface_it->get_adjacent_edges();
  size_t cnt = 0;
  for(vector<edgeit_type>::const_iterator adj_edge_it = adjacent_edges.begin();
      adj_edge_it != adjacent_edges.end(); ++adj_edge_it, ++cnt)
    {
      uv_values(zjucad::matrix::colon(), cnt) = ((*adj_edge_it)->get_uv_property()).get_uv();
    }
  return 0;
}

int Subdivision::get_edges_subdivided(const trifaceit_type &triface_it,
                                      vector<bool> &edge_flag)
{
  assert(edge_flag.size() == 3);
  int cnt = 0;
  const vector<edgeit_type> &adjacent_edges = triface_it->get_adjacent_edges();
  for(size_t i = 0; i < adjacent_edges.size(); ++i)
    {
      if(adjacent_edges[i]->get_size_field() > threshhold_)
      {
        edge_flag[i] = true;
        ++cnt;
      }
      else
        edge_flag[i] =false;
    }
  return cnt;
}

int Subdivision::sub_one_edge(const trifaceit_type &triface_it, const vector<bool> &edge_flag)
{
  adjust_points_index(triface_it, edge_flag);
  zjucad::matrix::matrix<double> face_frames(3, 3);
  zjucad::matrix::matrix<double> uv_values(2, 3);
  if(evfile_flag_)
   {

      get_face_three_frames(triface_it, face_frames);
      get_face_uv_values(triface_it, uv_values);
   }
  const vector<size_t> &face = triface_it->get_face();
  const vector<edgeit_type> &adjacent_edges = triface_it->get_adjacent_edges();
  double a[9];
  a[0] = 0.0;  a[1] = 0.5;  a[2] = 0.5;
  a[3] = 0.25; a[4] = 0.0;  a[5] = 0.75;
  a[6] = 0.25; a[7] = 0.5;  a[8] = 0.25;
  insert_new_face(triface_it, face[0], face[1], adjacent_edges[1]->get_mid_index(), 1, a, uv_values, face_frames);

  a[0] = 0.25; a[1] = 0.5; a[2] = 0.25;
  a[3] = 0.75; a[4] = 0.0; a[5] = 0.25;
  a[6] = 0.5 ; a[7] = 0.5; a[8] = 0.0;
  insert_new_face(triface_it, face[0], adjacent_edges[1]->get_mid_index(), face[2], 1, a, uv_values, face_frames);
  return 0;
}

int Subdivision::sub_two_edge(const trifaceit_type &triface_it, const vector<bool> &edge_flag)
{
  adjust_points_index(triface_it, edge_flag);
  zjucad::matrix::matrix<double> face_frames(3, 3);
  zjucad::matrix::matrix<double> uv_values(2, 3);
  if(evfile_flag_)
   {
      get_face_three_frames(triface_it, face_frames);
      get_face_uv_values(triface_it, uv_values);
   }
  const vector<size_t> &face = triface_it->get_face();
  const vector<edgeit_type> &adjacent_edges = triface_it->get_adjacent_edges();
  double a[9];
  a[0] = 0.0;  a[1] = 0.75; a[2] = 0.25;
  a[3] = 0.25; a[4] = 0.5;  a[5] = 0.25;
  a[6] = 0.25; a[7] = 0.75; a[8] = 0.0;
  insert_new_face(triface_it, face[0], adjacent_edges[0]->get_mid_index(), adjacent_edges[2]->get_mid_index(),
                  2, a, uv_values, face_frames);

  if(zjucad::matrix::norm(vertex_vec_[face[1]].get_coord() - vertex_vec_[adjacent_edges[2]->get_mid_index()].get_coord()) <
     zjucad::matrix::norm(vertex_vec_[face[2]].get_coord() - vertex_vec_[adjacent_edges[0]->get_mid_index()].get_coord()) )
    {
      a[0] = 0.0;  a[1] = 0.25; a[2] = 0.75;
      a[3] = 0.25; a[4] = 0.25; a[5] = 0.5;
      a[6] = 0.25; a[7] = 0.5;  a[8] = 0.25;
      insert_new_face(triface_it, adjacent_edges[0]->get_mid_index(), face[1], adjacent_edges[2]->get_mid_index(),
                      2, a, uv_values, face_frames);

      a[0] = 0.25; a[1] = 0.25; a[2] = 0.5;
      a[3] = 0.5;  a[4] = 0.0;  a[5] = 0.5;
      a[6] = 0.75; a[7] = 0.25; a[8] = 0.0;
      insert_new_face(triface_it, adjacent_edges[2]->get_mid_index(), face[1], face[2],
                      2, a, uv_values, face_frames);
    }
  else
    {
      a[0] = 0.0;  a[1] = 0.25; a[2] = 0.75;
      a[3] = 0.5;  a[4] = 0.0;  a[5] = 0.5;
      a[6] = 0.5;  a[7] = 0.25; a[8] = 0.25;
      insert_new_face(triface_it, adjacent_edges[0]->get_mid_index(), face[1], face[2],
                      2, a, uv_values, face_frames);

      a[0] = 0.5;  a[1] = 0.25; a[2] = 0.25;
      a[3] = 0.75; a[4] = 0.25; a[5] = 0.0;
      a[6] = 0.25; a[7] = 0.5;  a[8] = 0.25;
      insert_new_face(triface_it, adjacent_edges[0]->get_mid_index(), face[2], adjacent_edges[2]->get_mid_index(),
                      2, a, uv_values, face_frames);
    }

  return 0;
}

int Subdivision::sub_three_edge(const trifaceit_type &triface_it, const vector<bool> &edge_flag)
{
  adjust_points_index(triface_it, edge_flag);
  zjucad::matrix::matrix<double> face_frames(3, 3);
  zjucad::matrix::matrix<double> uv_values(2, 3);
  if(evfile_flag_)
   {
      get_face_three_frames(triface_it, face_frames);
      get_face_uv_values(triface_it, uv_values);
   }
  const vector<size_t> &face = triface_it->get_face();
  const vector<edgeit_type> &adjacent_edges = triface_it->get_adjacent_edges();
  double a[9];
  a[0] = 0.0;  a[1] = 0.75; a[2] = 0.25;
  a[3] = 0.25; a[4] = 0.5;  a[5] = 0.25;
  a[6] = 0.25; a[7] = 0.75; a[8] = 0.0;
  insert_new_face(triface_it, face[0], adjacent_edges[0]->get_mid_index(), adjacent_edges[2]->get_mid_index(),
                  3, a, uv_values, face_frames);

  a[0] = 0.0;  a[1] = 0.25; a[2] = 0.75;
  a[3] = 0.25; a[4] = 0.0;  a[5] = 0.75;
  a[6] = 0.25; a[7] = 0.25; a[8] = 0.5;
  insert_new_face(triface_it, adjacent_edges[0]->get_mid_index(), face[1], adjacent_edges[1]->get_mid_index(),
                  3, a, uv_values, face_frames);

  a[0] = 0.25; a[1] = 0.25; a[2] = 0.5;
  a[3] = 0.5;  a[4] = 0.25; a[5] = 0.25;
  a[6] = 0.25; a[7] = 0.5;  a[8] = 0.25;
  insert_new_face(triface_it, adjacent_edges[0]->get_mid_index(), adjacent_edges[1]->get_mid_index(), adjacent_edges[2]->get_mid_index(),
                  3, a, uv_values, face_frames);

  a[0] = 0.5;  a[1] = 0.25; a[2] = 0.25;
  a[3] = 0.75; a[4] = 0.0;  a[5] = 0.25;
  a[6] = 0.75; a[7] = 0.25; a[8] = 0.0;
  insert_new_face(triface_it, adjacent_edges[2]->get_mid_index(), adjacent_edges[1]->get_mid_index(), face[2],
                  3, a, uv_values, face_frames);
  return 0;
}

int Subdivision::adjust_points_index(const trifaceit_type &triface_it, const vector<bool> &edge_flag)
{
  vector<edgeit_type> &adjacent_edges = triface_it->get_adjacent_edges();
  vector<size_t> &face = triface_it->get_face();
  for(size_t i = 0; i < adjacent_edges.size(); ++i)
    if(edge_flag[i])
      insert_mid_point(adjacent_edges[i]);
  if(edge_flag[0] && !edge_flag[1] && !edge_flag[2]) // for one edge subdivided, the edge whose index is 1 is the one we want
    {
      swap(face[0], face[2]);
      swap(face[1], face[2]);
      swap(adjacent_edges[0], adjacent_edges[2]);
      swap(adjacent_edges[1], adjacent_edges[2]);
    }
  else if(!edge_flag[0] && !edge_flag[1] && edge_flag[2])
    {
      swap(face[0], face[1]);
      swap(face[1], face[2]);
      swap(adjacent_edges[0], adjacent_edges[1]);
      swap(adjacent_edges[1], adjacent_edges[2]);
    }

  else if(edge_flag[0] && edge_flag[1] && !edge_flag[2]) // for two edge subdivided, the edges whose index is 0 and 2 is we want
    {
      swap(face[0], face[1]);
      swap(face[1], face[2]);
      swap(adjacent_edges[0], adjacent_edges[1]);
      swap(adjacent_edges[1], adjacent_edges[2]);
    }
  else if(!edge_flag[0] && edge_flag[1] && edge_flag[2])
    {
      swap(face[0], face[2]);
      swap(face[1], face[2]);
      swap(adjacent_edges[0], adjacent_edges[2]);
      swap(adjacent_edges[1], adjacent_edges[2]);
    }
  return 0;
}

int Subdivision::update_adjacent_edges(const vector<size_t> &face,
                                       vector<edgeit_type> &adjacent_edges)
{
  assert(face.size() == 3 && adjacent_edges.size() == 3);
  pair<size_t, size_t> edge_pair;
  edge2it_type::iterator edge2it_it;
  for(size_t i = 0; i < face.size(); ++i)
    {
      edge_pair = make_pair(face[i], face[(i+1) % 3]);
      if(edge_pair.first > edge_pair.second)
        swap(edge_pair.first, edge_pair.second);
      edge2it_it = edge2it_.find(edge_pair);
      assert(edge2it_it != edge2it_.end());
      adjacent_edges[i] = edge2it_it->second;
      assert(edge_pair == adjacent_edges[i]->get_edge());
    }
  return 0;
}

int Subdivision::insert_mid_point(const edgeit_type &edge_it)
{
  if(edge_it->get_mid_index() == 0)
   {
      const pair<size_t, size_t> &edge_pair = edge_it->get_edge();
      zjucad::matrix::matrix<double> mid_point = (vertex_vec_[edge_pair.first].get_coord() +
                                                  vertex_vec_[edge_pair.second].get_coord()) / 2;
      vertex_vec_.push_back(Vertex(mid_point));
      edge_it->set_mid_index(vertex_vec_.size() - 1);
    }
  return 0;
}

int Subdivision::insert_new_face(const trifaceit_type &triface_it,
                                 const size_t id0, const size_t id1, const size_t id2,
                                 const size_t sub_edge_num, const double *barycentric,
                                 const zjucad::matrix::matrix<double> &uv_value,
                                 const zjucad::matrix::matrix<double> &frame_value)
{
  triface_list_.push_back(Tri_face(id0, id1, id2));
  trifaceit_type cur_face_it = --triface_list_.end();
  vector<edgeit_type> &adjacent_edges = cur_face_it->get_adjacent_edges();
  assert(adjacent_edges.size() == 0);
  edgeit_type cur_edge_it;
  vector<size_t> temp_vec(3, 0);
  temp_vec[0] = id0;
  temp_vec[1] = id1;
  temp_vec[2] = id2;
  const vector<size_t> &face = triface_it->get_face();
  vector<bool> is_in(3, false);
  for(size_t i = 0; i < temp_vec.size(); ++i)
    if(is_element_in_vec(temp_vec[i], face))
      is_in[i] = true;
  bool flag;
  for(size_t i = 0; i < 3; ++i)
    {
      if(is_in[i] && is_in[(i+1)%3])
        flag =false;
      else
        flag =true;
      insert_new_edge(temp_vec[i], temp_vec[(i+1)%3], flag, cur_face_it, triface_it,
                      uv_value, frame_value, barycentric + 3*i, cur_edge_it);
      adjacent_edges.push_back(cur_edge_it);
    }
  assert(adjacent_edges.size() == 3);
  return 0;
}

int Subdivision::insert_new_edge(const size_t id0, const size_t id1, const bool flag,
                                 const trifaceit_type &cur_face_it,
                                 const trifaceit_type &old_face_it,
                                 const zjucad::matrix::matrix<double> &uv_value,
                                 const zjucad::matrix::matrix<double> &frame_value,
                                 const double *barycentric, edgeit_type &cur_edge_it)
{
  assert(id0 != id1);
  pair<size_t, size_t> edge_pair(id0, id1);
  if(id0 > id1)
    swap(edge_pair.first, edge_pair.second);
  const edge2it_type::iterator &edge2it_it = edge2it_.find(edge_pair);
  if(edge2it_it == edge2it_.end())
    {
      edge_list_.push_back(Edge(id0, id1));
      cur_edge_it = --edge_list_.end();
      edge2it_.insert(make_pair(edge_pair, cur_edge_it));
      cur_edge_it->set_length(zjucad::matrix::norm(vertex_vec_[id0].get_coord() -
                                                   vertex_vec_[id1].get_coord()));
      edgeframe_list_type &edgeframe_list = cur_edge_it->get_edgeframe_list();
      assert(edgeframe_list.size() == 0);
      zjucad::matrix::matrix<double> weight(3, 1);
      weight[0] = barycentric[1] + barycentric[2] - barycentric[0];
      weight[1] = barycentric[0] + barycentric[2] - barycentric[1];
      weight[2] = barycentric[0] + barycentric[1] - barycentric[2];
      edgeframe_list.push_back(make_pair(cur_face_it, Frame_property()));
      edgeframe_list_type::iterator edgeframe_it = --edgeframe_list.end();
      if(evfile_flag_)
        {
          (edgeframe_it->second).cal_property(weight, frame_value);
          cur_edge_it->get_uv_property().cal_property(weight, uv_value);
        }
    }
  else
    {
      cur_edge_it = edge2it_it->second;
      if(flag)
        {
          edgeframe_list_type &edgeframe_list = cur_edge_it->get_edgeframe_list();
          assert(edgeframe_list.size() != 0);
          edgeframe_list.push_back(make_pair(cur_face_it, Frame_property()));
          edgeframe_list_type::iterator edgeframe_it = --edgeframe_list.end();
          zjucad::matrix::matrix<double> weight(3, 1);
          weight[0] = barycentric[1] + barycentric[2] - barycentric[0];
          weight[1] = barycentric[0] + barycentric[2] - barycentric[1];
          weight[2] = barycentric[0] + barycentric[1] - barycentric[2];
          if(evfile_flag_)
            (edgeframe_it->second).cal_property(weight, frame_value);
        }
      else
        {
          edgeframe_list_type &edgeframe_list = cur_edge_it->get_edgeframe_list();
          assert(edgeframe_list.size() != 0);
          for(edgeframe_list_type::iterator edgeframe_it = edgeframe_list.begin();
              edgeframe_it != edgeframe_list.end(); ++edgeframe_it)
            {
              if(edgeframe_it->first->get_face() == old_face_it->get_face())
                {
                  edgeframe_it->first = cur_face_it;
                }
            }
        }
    }
  return 0;
}

int Subdivision::load_flfile(const char *fl_file)
{
  if(!fl_file)
    return 2;
  std::ifstream ifs(fl_file);
  if(!ifs)
    {
      std::cerr<<"load FL file error!"<<std::endl;
      return 1;
    }
  size_t num, len, temp;
  ifs >> num;
  feature_line_.resize(num);
  for(size_t i = 0; i < num; ++i)
    {
      ifs >> len;
      for(size_t j = 0; j < len; ++j)
        {
          ifs >> temp;
          feature_line_[i].push_back(temp);
        }
    }
  ifs.close();
  std::cout<< "load FL file end." << std::endl;
  return 0;
}

int Subdivision::save_flfile(const char *fl_file)
{
  std::ofstream ofs(fl_file);
  if( !ofs )
    return 1;
  size_t num=feature_line_.size();
  ofs << num << endl;
  for(size_t i = 0; i < num; ++i)
    {
      ofs << feature_line_[i].size() << endl;
      for(list<size_t>::const_iterator it = feature_line_[i].begin();
          it != feature_line_[i].end(); ++it)
        {
          ofs << *it << " ";
        }
      ofs << endl;
    }
  ofs.close();
  cout << "save FL file end." << endl;
  return 0;
}

int Subdivision::update_feature_line()
{
  pair<size_t, size_t> edge_pair;
  edge2it_type::iterator edge2it_it;
  size_t mid_index(0);
  for(size_t i = 0; i < feature_line_.size(); ++i)
    {
      for(list<size_t>::iterator it = feature_line_[i].begin(), it1; ;)
        {
          it1 = it;
          ++it1;
          if(it1 == feature_line_[i].end())
            break;
          edge_pair = make_pair(*it, *it1);
          if(edge_pair.first > edge_pair.second)
            swap(edge_pair.first, edge_pair.second);
          edge2it_it = edge2it_.find(edge_pair);
          assert(edge2it_it != edge2it_.end());
          mid_index = edge2it_it->second->get_mid_index();
          if(mid_index == 0) // the egde was not subdivided
            {
              ++it;
            }
          else
            {
              feature_line_[i].insert(it1, mid_index);
            }
        }
    }
  return 0;
}

int Subdivision::delete_remain_edges()
{
  edgeit_type edge_it1, edge_it2;
  for(edge_it1 = edge_list_.begin(); edge_it1 != edge_list_.end();)
    {
      const edgeframe_list_type &edgeframe_list = edge_it1->get_edgeframe_list();
      if(edgeframe_list.size() == 0)
        {
          edge_it2 = edge_it1++;
          edge_list_.erase(edge_it2);
        }
      else
        ++edge_it1;
    }
  return 0;
}

int Subdivision::import_mesh(const zjucad::matrix::matrix<double> &points,
                             const zjucad::matrix::matrix<size_t> &tri_face)
{
  assert(points.size(1) == 3 && tri_face.size(1) == 3);
  for(size_t i = 0; i < points.size(2); ++i)
    {
      vertex_vec_.push_back(Vertex(points(zjucad::matrix::colon(), i)));
    }
  trifaceit_type cur_face_it;
  for(size_t i = 0; i < tri_face.size(2); ++i)
    {
      add_face(tri_face(0, i), tri_face(1, i), tri_face(2, i), cur_face_it);
      cur_face_it->set_face_index(i);
    }
  esfile_flag_ = false;
  evfile_flag_ = false;
  flfile_flag_ = false;
  return 0;
}

int Subdivision::export_mesh(zjucad::matrix::matrix<double> &points,
                             zjucad::matrix::matrix<size_t> &tri_face)
{
  points.resize(3, vertex_vec_.size());
  tri_face.resize(3, triface_list_.size());
  for(size_t i = 0; i < vertex_vec_.size(); ++i)
    {
      points(zjucad::matrix::colon(), i) = vertex_vec_[i].get_coord();
    }

  size_t i = 0;
  for(trifaceit_type triface_it = triface_list_.begin();
      triface_it != triface_list_.end(); ++triface_it, ++i)
    {
      const vector<size_t> &face = triface_it->get_face();
      for(size_t j = 0; j < 3; ++j)
        {
          tri_face(j, i) = face[j];
        }
    }
  return 0;
}

int Subdivision::subdived_edges(const vector<pair<size_t, size_t> > &sub_edges, vector<size_t> &insert_points)
{
  pair<size_t, size_t> edge_pair;
  edge2it_type::iterator edge2it_it;
  const double flag1 = 100, flag2 = 50;
  for(size_t i = 0; i < sub_edges.size(); ++i)
    {
      edge_pair = sub_edges[i];
      if(edge_pair.first > edge_pair.second)
        swap(edge_pair.first, edge_pair.second);
      edge2it_it = edge2it_.find(edge_pair);
      assert(edge2it_it != edge2it_.end());
      edge2it_it->second->set_size_field(flag1);
    }
  set_threshhold(flag2);
  cout << "threshhold:" << threshhold_ << endl;
  flfile_flag_ = true;   //do not delete edges
  adaptive_subdivided();
  insert_points.resize(sub_edges.size());
  for(size_t i = 0; i < sub_edges.size(); ++i)
    {
      edge_pair = sub_edges[i];
      if(edge_pair.first > edge_pair.second)
        swap(edge_pair.first, edge_pair.second);
      edge2it_it = edge2it_.find(edge_pair);
      assert(edge2it_it != edge2it_.end());
      insert_points[i] = edge2it_it->second->get_mid_index();
    }
  return 0;
}

int Subdivision::flip_edges(const vector<pair<size_t, size_t> > &edges)
{
  pair<size_t, size_t> edge_pair;
  edge2it_type::iterator edge2it_it;
  flfile_flag_ = false;
  for(size_t i = 0; i < edges.size(); ++i)
    {
      edge_pair = edges[i];
      if(edge_pair.first > edge_pair.second)
        swap(edge_pair.first, edge_pair.second);
      edge2it_it = edge2it_.find(edge_pair);
      assert(edge2it_it != edge2it_.end());
      edge_flip(edge2it_it->second);
    }
  return 0;
}


int Subdivision::edge_flip(const edgeit_type &edge_it)
{
  assert(edge_it != edge_list_.end());
  size_t len = edge_it->get_edgeframe_list().size();
  if(len != 2)
    {
      cerr << "the edge can not do this operation" << endl;
      return 1;
    }
  pair<size_t, size_t> edge_pair = edge_it->get_edge();
  vector<size_t> new_points(2, 0);
  vector<trifaceit_type> del_face_vec;
  const edgeframe_list_type &edgeframe_list = edge_it->get_edgeframe_list();
  for(edgeframe_list_type::const_iterator edgeframe_it = edgeframe_list.begin();
      edgeframe_it != edgeframe_list.end(); ++edgeframe_it)
    {
      del_face_vec.push_back(edgeframe_it->first);
    }
  for(size_t i = 0; i < del_face_vec.size(); ++i)
    {
      const vector<size_t> &face = del_face_vec[i]->get_face();
      new_points[i] = face[0] + face[1] + face[2] - edge_pair.first - edge_pair.second;
      delete_face(del_face_vec[i]);
    }
  trifaceit_type triface_it;
  add_face(new_points[0], edge_pair.first, new_points[1], triface_it);
  add_face(new_points[0], new_points[1], edge_pair.second, triface_it);
  return 0;
}
