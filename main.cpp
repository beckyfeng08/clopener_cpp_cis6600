#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <chrono>
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
          : std::string("../data/godzilla.obj");

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
  std::unique_ptr<ClosingFlow> flow;
  bool running = false;
  std::chrono::steady_clock::time_point flow_start_time;

  // ~~~~~~~~~~~ STEP 2: HANDLE INPUTS ~~~~~~~~~~~~

  std::cerr << "Press c to run closing_flow on the current mesh.\n";

  viewer.callback_key_pressed =
      [&](igl::opengl::glfw::Viewer &v, unsigned int key, int) -> bool {
        if (key != 'c' && key != 'C') return false;

        if (running) {
          std::cerr << "Flow already running, ignoring key press\n";
          return true;
        }

        igl::write_triangle_mesh("../data/iters/iter_0000.obj", V, F);

        try {
          flow = std::make_unique<ClosingFlow>(V, F, params);
          running = true;
          flow_start_time = std::chrono::steady_clock::now();   // start timer
          v.core().is_animating = true;
          std::cerr << "Started closing_flow\n";
        } catch (const std::exception& e) {
          std::cerr << "Failed to start closing_flow: " << e.what() << "\n";
        }
        return true;
      };

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

        // Stop condition
        if (!keep_going || flow->iteration() >= params.maxiter) {
          running = false;
          v.core().is_animating = false;

          // ~~~~~~~~~~~ STEP 4: STOP TIMER & REPORT ~~~~~~~~~~~~
          const double total_seconds =
              std::chrono::duration<double>(
                  std::chrono::steady_clock::now() - flow_start_time).count();
          const double remesh_seconds  = flow->remesh_seconds_total();
          const double other_seconds   = total_seconds - remesh_seconds;
          const int    n_iters         = flow->iteration();
          const double per_iter        = (n_iters > 0)
                                       ? total_seconds / n_iters
                                       : 0.0;

          if (!igl::write_triangle_mesh("../data/mesh_remeshed.obj", V, F)) {
            std::cerr << "Failed to write ../data/mesh_remeshed.obj\n";
          } else {
            std::cerr << "Successfully wrote ../data/mesh_remeshed.obj\n";
          }
          std::cerr << "closing_flow finished after " << n_iters
                    << " iterations (snapshots saved as iter_NNNN.obj)\n";
          std::cerr << "  total time:    " << total_seconds  << " s\n";
          std::cerr << "  remesh time:   " << remesh_seconds << " s\n";
          std::cerr << "  other time:    " << other_seconds  << " s\n";
          std::cerr << "  avg per iter:  " << per_iter       << " s\n";
          std::cerr << "Press c to run again on current mesh.\n";
        }

        return false;
      };

  viewer.launch();
  return 0;
}