#include <cstdio>
#include <cassert>
#include <cstring>

#include "basic_define.h"

//size_t sxx::hash_value(const Edge &edge)
//{
//  boost::hash<std::pair<size_t, size_t> > t_hash;
//  return t_hash(edge.get_edge());
//}

//size_t sxx::hash_value(const Tri_face &face)
//{
//  boost::hash<std::vector<size_t> > t_hash;
//  return t_hash(face.get_face());
//}

int sxx::load_obj(const char * file_name,
                  std::vector<std::vector<double> > &points,
                  std::vector<size_t> &faces)
{
  const size_t BUF_SIZE = 128;
  FILE *file = fopen(file_name,"r");
  if(!file)
    return 1;
  char buf[BUF_SIZE];
  std::vector<double> ver(3, 0.0), trash;
  while (fscanf(file,"%s",buf) != EOF)
    {
      switch (buf[0])
        {
        case '#': /* comment */
          // eat up rest of line
          fgets(buf,sizeof(buf),file);
          break;
        case 'v': // v, vn, vt
          switch (buf[1])
            {
            case '\0':			/* vertex */
              fscanf(file, "%lf %lf %lf", &(ver[0]), &(ver[1]), &(ver[2]));
              points.push_back(ver);
              break;
            case 'n':				/* normal */
              fscanf(file, "%lf %lf %lf", &(trash), &(trash), &(trash));
              break;
            case 't':				/* texcoord */
              for (int j = 0;j < 2; ++j)   //only accept 2 component texture coordinate
                fscanf(file,"%lf",&(trash));
              break;
            default:
              assert(false);
              goto errorEnd;
            }
          break;
        case 'm':  //ignore material
          fgets(buf, sizeof(buf), file);
          sscanf(buf, "%s %s", buf, buf);
          //model->mtllibname = strdup(buf);
          //ReadMaterial(buf,mesh);
          break;
        case 'u':  //use material, ignore
          /* eat up rest of line */
          fgets(buf, sizeof(buf), file);
          break;
        case 'g':  //group, ignore
          /* eat up rest of line */
          fgets(buf, sizeof(buf), file);
          sscanf(buf, "%s", buf);
          break;
        case 's': // Smoothing group.ignore
          fgets(buf,sizeof(buf),file);
          sscanf(buf,"%s",buf);
          break;
        case 'f':				// face
          {
            int v = 0, n = 0, t = 0;
            fscanf(file, "%s", buf);
            /* can be one of %d, %d//%d, %d/%d, %d/%d/%d %d//%d */
            if (strstr(buf, "//"))
              {
                /* v//n */
                sscanf(buf, "%d//%d", &v, &n); // Indices start with 1 in obj, so:
                if(v < 0)
                  v = points.size() + v + 1;
                faces.push_back(v-1);

                fscanf(file, "%d//%d", &v, &n);
                if(v < 0)
                  v = points.size() + v + 1;
                faces.push_back(v-1);

                fscanf(file, "%d//%d", &v, &n);
                if(v < 0)
                  v = points.size() + v + 1;
                faces.push_back(v-1);
              }
            else if (sscanf(buf, "%d/%d/%d", &v, &t, &n) == 3)
              {
                /* v/t/n */
                if(v < 0)
                  v = points.size() + v + 1;
                faces.push_back(v-1);

                fscanf(file, "%d/%d/%d", &v, &t, &n);
                if(v < 0)
                  v = points.size() + v + 1;
                faces.push_back(v-1);

                if(v < 0)
                  v = points.size() + v + 1;
                  faces.push_back(v-1);
              }
            else if (sscanf(buf, "%d/%d", &v, &t) == 2)
              {
                /* v/t */
                if(v < 0)
                  v = points.size() + v + 1;
                faces.push_back(v-1);

                fscanf(file, "%d/%d", &v, &t);
                if(v < 0)
                  v = points.size() + v + 1;
                faces.push_back(v-1);

                fscanf(file, "%d/%d", &v, &t);
                if(v < 0)
                  v = points.size() + v + 1;
                faces.push_back(v-1);
              }
            else
              {
                /* v */
                sscanf(buf, "%d", &v);
                if(v < 0)
                  v = points.size() + v + 1;
                faces.push_back(v-1);

                fscanf(file,"%d", &v);
                if(v < 0)
                  v = points.size() + v + 1;
                faces.push_back(v-1);

                fscanf(file,"%d", &v);
                if(v < 0)
                  v = points.size() + v + 1;
                faces.push_back(v-1);
              }
          }
          break;
        default:
          // eat up rest of line
          fgets(buf,sizeof(buf),file);
        }
    }

errorEnd:
  fclose(file);
  return 0 ;

}
