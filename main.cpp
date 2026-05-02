#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <cstdio>
#include <iostream>
#include <memory>
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
  } else {
    std::cerr << "Successfully wrote ../data/mesh_initial.obj\n";
  }

  ClosingFlowParams params;

  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);

  // ~~~~~~~~~~~ STEP 1: SHARED STATE FOR THE STEPPER ~~~~~~~~~~~~
  // These are owned by main() and captured by reference in the callbacks.
  std::unique_ptr<ClosingFlow> flow;
  bool running = false;

  // ~~~~~~~~~~~ STEP 2: HANDLE INPUTS ~~~~~~~~~~~~

  std::cerr << "Press c to run closing_flow on the current mesh.\n";

  // 'c' key starts a new flow. The actual stepping happens in the pre-draw callback.
  viewer.callback_key_pressed =
      [&](igl::opengl::glfw::Viewer &v, unsigned int key, int) -> bool {
        if (key != 'c' && key != 'C') {
          return false;
        }

        if (running) {
          std::cerr << "Flow already running, ignoring key press\n";
          return true;
        }

        // Snapshot iteration 0 (initial state) so the user can scrub back to it
        igl::write_triangle_mesh("../data/iters/iter_0000.obj", V, F);

        try {
          flow = std::make_unique<ClosingFlow>(V, F, params);
          running = true;
          v.core().is_animating = true;  // make the viewer redraw continuously
          std::cerr << "Started closing_flow\n";
        } catch (const std::exception& e) {
          std::cerr << "Failed to start closing_flow: " << e.what() << "\n";
        }
        return true;
      };

  // pre_draw fires every frame the viewer renders. While `running` is true,
  // we run one closing_flow iteration per frame, push the new mesh to the
  // viewer, and write the snapshot to disk.
  viewer.callback_pre_draw =
      [&](igl::opengl::glfw::Viewer &v) -> bool {
        if (!running || !flow) return false;

        bool keep_going = false;
        try {
          keep_going = flow->step();
        } catch (const std::exception& e) {
          std::cerr << "closing_flow step failed: " << e.what() << "\n";
          running = false;
          v.core().is_animating = false;
          return false;
        }

        // Pull the updated mesh and push to the viewer
        V = flow->current_V();
        F = flow->current_F();
        v.data().clear();
        v.data().set_mesh(V, F);
        v.data().set_face_based(true);

        // ~~~~~~~~~~~ STEP 3: SNAPSHOT THIS ITERATION TO DISK ~~~~~~~~~~~~
        char fname[256];
        std::snprintf(fname, sizeof(fname),
                      "../data/iters/iter_%04d.obj", flow->iteration());
        if (!igl::write_triangle_mesh(fname, V, F)) {
          std::cerr << "Failed to write " << fname << "\n";
        }

        // Stop condition: flow says we're done, or we've hit maxiter
        if (!keep_going || flow->iteration() >= params.maxiter) {
          running = false;
          v.core().is_animating = false;

          if (!igl::write_triangle_mesh("../data/mesh_remeshed.obj", V, F)) {
            std::cerr << "Failed to write ../data/mesh_remeshed.obj\n";
          } else {
            std::cerr << "Successfully wrote ../data/mesh_remeshed.obj\n";
          }
          std::cerr << "closing_flow finished after " << flow->iteration()
                    << " iterations (snapshots saved as iter_NNNN.obj)\n";
          std::cerr << "Press c to run again on current mesh.\n";
        }

        return false;  // don't suppress the actual draw
      };

  viewer.launch();
  return 0;
}