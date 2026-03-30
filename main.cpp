#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <iostream>
#include <string>
#include "closing_flow.h"

int main(int argc, char *argv[])
{
  // ~~~~~~~~~~~ STEP 0: LOAD MESH ~~~~~~~~~~~~
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  const std::string mesh_path =
      (argc > 1 && argv[1] != nullptr && argv[1][0] != '\0')
          ? std::string(argv[1])
          : std::string("../data/bunny.obj");

  if (!igl::read_triangle_mesh(mesh_path, V, F)) {
    std::cerr << "Failed to load mesh: " << mesh_path << "\n";
    std::cerr << "Usage: " << (argc > 0 ? argv[0] : "example")
              << " [mesh_path]\n";
    return 1;
  }

  if (!igl::write_triangle_mesh("../data/mesh_initial.obj", V, F)) {
    std::cerr << "Failed to write ../data/mesh_initial.obj\n";
  } else{
    std::cerr << "Successfully wrote ../data/mesh_initial.obj\n";
  }

  ClosingFlowParams params;

  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);

  // ~~~~~~~~~~~ STEP 1: HANDLE INPUTS ~~~~~~~~~~~~

  std::cerr << "Press c to run closing_flow on the current mesh.\n";

  viewer.callback_key_pressed =
      [&](igl::opengl::glfw::Viewer &v, unsigned int key, int) -> bool {
        if (key != 'c' && key != 'C') {
          return false;
        }

        // ~~~~~~~~~~~ STEP 2: CLOSE THE MESH ~~~~~~~~~~~~

        Eigen::MatrixXd Vout;
        Eigen::MatrixXi Fout;

        bool closed = closing_flow(V, F, params, Vout, Fout); // running closing operations on original vertices and faces

        if (!closed) {
          std::cerr << "closing_flow failed\n";
          return true;
        }
        V = std::move(Vout);
        F = std::move(Fout);
        // Remeshing changes |V|/|F|; ViewerData::set_mesh requires clear first
        v.data().clear();
        v.data().set_mesh(V, F);
        v.data().set_face_based(true);

          // ~~~~~~~~~~~ STEP 3: WRITE NEW MESH TO OUTPUT OBJ ~~~~~~~~~~~~

        if (!igl::write_triangle_mesh("../data/mesh_remeshed.obj", V, F)) {
          std::cerr << "Failed to write ../data/mesh_remeshed.obj\n";
        } else{
          std::cerr << "Successfully wrote ../data/mesh_remeshed.obj\n";
        }

        std::cerr << "closing_flow finished; mesh updated (press c to run again on "
                     "current mesh)\n";
        return true;
      };

  viewer.launch();
  return 0;
}
