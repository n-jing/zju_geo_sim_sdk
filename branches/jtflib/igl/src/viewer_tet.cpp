#include "../include/viewer_tet.h"
#include "../include/readtet.h"
#include <jtflib/mesh/mesh.h>

namespace igl{

  template <typename v1_type, typename v2_type>
  void eigen_matrix2zju_matrix(const Eigen::PlainObjectBase<v1_type> & A,
                               zjucad::matrix::matrix<v2_type> &B)
  {
    B.resize(A.cols(), A.rows());
    for(size_t ri = 0; ri < B.size(1); ++ri){
        for(size_t ci = 0; ci < B.size(2); ++ci){
            B(ri,ci) = A(ci,ri);
          }
      }
  }

  template <typename v1_type, typename v2_type>
  void zju_matrix2eigen_matrix(const zjucad::matrix::matrix<v1_type> &A,
                               Eigen::PlainObjectBase<v2_type> &B)
  {
    B.resize(A.size(2), A.size(1));
    for(size_t ri = 0; ri < B.rows(); ++ri){
        for(size_t ci = 0; ci < B.cols(); ++ci){
            B(ri,ci) = A(ci,ri);
          }
      }
  }

  bool Viewer_tet::load_tet_from_file(const char *mesh_file_name)
  {
    std::string mesh_file_name_string = std::string(mesh_file_name);

    // first try to load it with a plugin
    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->load(mesh_file_name_string))
        return true;

    data.clear();

    size_t last_dot = mesh_file_name_string.rfind('.');
    if (last_dot == std::string::npos)
      {
        printf("Error: No file extension found in %s\n",mesh_file_name);
        return false;
      }

    if (!igl::readtet(mesh_file_name_string, data.V, T))
      return false;
    zjucad::matrix::matrix<double> node;
    zjucad::matrix::matrix<size_t> tet;
    eigen_matrix2zju_matrix(data.V,node);
    eigen_matrix2zju_matrix(T,tet);
    zjucad::matrix::matrix<size_t> face;
    std::unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
    if(!fa.get()) throw std::logic_error("can not build face2tet_adjacent.");

    jtf::mesh::get_outside_face(*fa, face, true, &node);
    jtf::mesh::save_obj("test_face.obj", face, node);
    zju_matrix2eigen_matrix(face, data.F);

    data.compute_normals();
    data.uniform_colors(Eigen::Vector3d(51.0/255.0,43.0/255.0,33.3/255.0),
                        Eigen::Vector3d(255.0/255.0,228.0/255.0,58.0/255.0),
                        Eigen::Vector3d(255.0/255.0,235.0/255.0,80.0/255.0));
    if (data.V_uv.rows() == 0)
      data.grid_texture();

    core.align_camera_center(data.V,data.F);

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->post_load())
        return true;

    return true;
  }

  bool Viewer_tet::save_tet_to_file(const char *mesh_file_name)
  {

  }

  int Viewer_tet::launch(std::string filename)
  {
    GLFWwindow* window;

    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
      return EXIT_FAILURE;

    glfwWindowHint(GLFW_SAMPLES, 16);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    #ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif
    window = glfwCreateWindow(1280, 800, "IGL Viewer", NULL, NULL);
    if (!window)
    {
      glfwTerminate();
      return EXIT_FAILURE;
    }

  glfwMakeContextCurrent(window);

#ifndef __APPLE__
  glewExperimental = true;
  GLenum err = glewInit();
  if (GLEW_OK != err)
  {
    /* Problem: glewInit failed, something is seriously wrong. */
    fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
  }
  fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
#endif

    #ifdef DEBUG
      int major, minor, rev;
      major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
      minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
      rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
      printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
      printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
      printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));
    #endif

    glfwSetInputMode(window,GLFW_CURSOR,GLFW_CURSOR_NORMAL);

    // Initialize AntTweakBar
    TwInit(TW_OPENGL_CORE, NULL);
    TwCopyStdStringToClientFunc(static_cast<TwCopyStdStringToClient>(::copy_str));


    // Initialize IGL viewer
    init();
    __viewer = this;

    // Register callbacks
    glfwSetKeyCallback(window, glfw_key_callback); glfwSetCursorPosCallback(window,glfw_mouse_move); glfwSetWindowSizeCallback(window,glfw_window_size); glfwSetMouseButtonCallback(window,glfw_mouse_press); glfwSetScrollCallback(window,glfw_mouse_scroll); glfwSetCharCallback(window, glfw_char_callback);

    // Handle retina displays (windows and mac)
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    int width_window, height_window;
    glfwGetWindowSize(window, &width_window, &height_window);

    highdpi = width/width_window;

    glfw_window_size(window,width_window,height_window);

    opengl.init();

    // Load the mesh passed as input
    if (filename.size() > 0)
      load_tet_from_file(filename.c_str());

    core.align_camera_center(data.V,data.F);

    // Rendering loop
    while (!glfwWindowShouldClose(window))
    {
      double tic = get_seconds();
      draw();

      glfwSwapBuffers(window);
      if(core.is_animating)
      {
        glfwPollEvents();
        // In microseconds
        double duration = 1000000.*(get_seconds()-tic);
        const double min_duration = 1000000./core.animation_max_fps;
        if(duration<min_duration)
        {
          std::this_thread::sleep_for(std::chrono::microseconds((int)(min_duration-duration)));
        }
      }
      else
      {
        glfwWaitEvents();
      }
    }

    opengl.free();
    core.shut();

    shutdown_plugins();

    glfwDestroyWindow(window);
    glfwTerminate();
    return EXIT_SUCCESS;
  }

}
