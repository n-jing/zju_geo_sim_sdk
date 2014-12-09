#include "dl_petsc_support.h"
#include <iostream>

static int get_function_point(const char* file_name,const char* fun_name, gpointer *point){

	GModule *module = g_module_open (file_name, G_MODULE_BIND_LAZY);
        // printf("[zsw_info] %s:%d get fun: %s from file: %s\n", __FILE__, __LINE__, fun_name, file_name);
	if (!module)
	{
		// printf("%s\n",g_module_error());
		return 1;
	}
	if (!g_module_symbol (module, fun_name, point))
	{
	    g_module_close (module);
		return 2;
	}
	if (point == NULL)
	{
		g_module_close(module);
		return 3;
	}
	return 0;
}

PETsc* ucreate_PETsc(const char * path)
{
  PETsc* (*create_PETsc)() = NULL;  
  std::string filename = path ? path : std::string(LIB_PREFIX) + "dl_petsc" + LIB_SUFFIX;
  if( get_function_point(filename.c_str(), "create_PETsc", (gpointer*)&create_PETsc))
  {
    std::cerr << "can't get_function_point create_PETsc from file"
              << filename
              << std::endl;
    return NULL;
  }
  return create_PETsc();
}

PETsc_CG *ucreate_PETsc_CG(const hj::sparse::csc<double, int32_t> &A, const char *pc_str, const char* path)
  {
    PETsc_CG* (*create_PETsc_CG)(const hj::sparse::csc<double, int32_t>&, const char *) = NULL;
    std::string filename = path ? path : std::string(LIB_PREFIX) + "dl_petsc" + LIB_SUFFIX;
    if (get_function_point(filename.c_str(), "create_PETsc_CG", (gpointer*)&create_PETsc_CG))
    {
      std::cerr << "can't get_function_point create_PETsc_CG from file:" << filename << std::endl;
      return NULL;
    }
    return create_PETsc_CG(A, pc_str);
  }
