#ifndef READ_GMESH_H
#define READ_GMESH_H

#include <fstream>
#include <boost/shared_ptr.hpp>
#include <gmesh/mesh/mesh.h>

template<class geometry_t, class connectivity_t>
int read_from_file(boost::shared_ptr<gmesh::mesh<geometry_t, connectivity_t> > mesh_model, const std::string& file_name) {
  int prefix_len = file_name.rfind('.');
  std::string suffix = file_name.substr(prefix_len + 1, file_name.length() - prefix_len - 1);
#ifdef _DEBUG
  std::cout << "[Debug Info] : " << "file name is " << file_name << std::endl;
  std::cout << "[Debug Info] : " << "prefix length is " << prefix_len << std::endl;
  std::cout << "[Debug Info] : " << "model file format is " << suffix << std::endl;
#endif
  if(suffix == "obj")
    read_obj(mesh_model, file_name);
  else {
    std::cerr << "[Error] : No support for file format " << suffix << std::endl;
    return 1;
  }
  return 0;
}

template<class geometry_t, class connectivity_t>
int read_obj(boost::shared_ptr<gmesh::mesh<geometry_t, connectivity_t> > mesh_model, const std::string& file_name) {
  std::ifstream fin(file_name.c_str(), std::ios_base::in);
  if(!fin) {
    std::cerr << "Cannot read file : " << file_name << std::endl;
    return 1;
  }

  std::string line, key_word;
  std::vector<vert_handle> vhandls;
  double x, y, z;

  while(getline(fin, line)) {
    if(fin.bad()) {
      std::cerr << "Warning! cannot read file properly!" << std::endl;
      return 1;
    }
    string_trim(line);
    if(line.length() == 0 || line[0] == '#') continue;

    std::stringstream stream(line);
    stream >> key_word;
    if(key_word == "v") {
      stream >> x >> y >> z;
      if(!stream.fail()) {
        mesh_model->add_vert(point(x, y, z));
      }
    }else if(key_word == "vt"){
      ///!TODO: vertex texture
    }else if(key_word == "vn"){
      ///!TODO: vertex normal
    }else if(key_word == "vc") {
      ///!TODO: vertex color
    }else if(key_word == "f") {
      int component(0), nv(0), index;
      vhandls.clear();
      std::string face_line;
      std::getline(stream, face_line);
      std::stringstream line_data(face_line);

      while(!line_data.eof()) {
        std::string vert_data;
        line_data >> vert_data;
        do {
          size_t found = vert_data.find("/");
          if(found != std::string::npos) {
            std::stringstream tmp(vert_data.substr(0, found));
            if(vert_data.substr(0, found).empty()) {
              vert_data = vert_data.substr(found+1);
              ++component;
              continue;
            }
            tmp >> index;
            vert_data = vert_data.substr(found+1);
          }else {
            std::stringstream tmp(vert_data);
            tmp >> index;
            vert_data = "";
            if(tmp.fail()) continue;
          }
          switch(component) {
            case 0: // vertex
              if(index < 0) {
                index = mesh_model->get_vert_num() + index + 1;
              }
              vhandls.push_back(vert_handle(index-1));
              break;
            case 1: // texture
              ///TODO:
              break;
            case 2: // normal
              ///TODO:
              break;
            default:
              break;
          }
          ++component;
        } while(!vert_data.empty());
        component=0;
        nv ++;
      }

      mesh_model->add_face(vhandls);
    }
  }
  fin.close();
  return 0;
}

#endif
