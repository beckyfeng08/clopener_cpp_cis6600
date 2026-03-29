#include <igl/opengl/glfw/Viewer.h>
#include "remesh_botsch.h"
#include <igl/boundary_loop.h>
#include "all_boundary_loop.h"
#include <Eigen/Core>
// #include <igl/read_triangle_mesh.h>
// #include <igl/write_triangle_mesh.h>


void create_cube(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
V= (Eigen::MatrixXd(8,3)<<
    0.0,0.0,0.0,
    0.0,0.0,1.0,
    0.0,1.0,0.0,
    0.0,1.0,1.0,
    1.0,0.0,0.0,
    1.0,0.0,1.0,
    1.0,1.0,0.0,
    1.0,1.0,1.0).finished();
  F = (Eigen::MatrixXi(12,3)<<
    0,6,4,
    0,2,6,
    0,3,2,
    0,1,3,
    2,7,6,
    2,3,7,
    4,6,7,
    4,7,5,
    0,4,5,
    0,5,1,
    1,5,7,
    1,7,3).finished();
}
int main(int argc, char *argv[])
{
  //  ~~~~~ STEP 0: LOAD MESH ~~~~~
  // Inline mesh of a cube

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  // generates a cube mesh and populates V and F
  create_cube(V, F);

  // uncomment to generate bunny mesh or any arbitrary mesh
  // igl::read_triangle_mesh(argc>1 ? argv[1] : "../data/bunny.off", V, F);

  // uncomment for read .obj files, be sure to uncomment headers too
      // igl::read_triangle_mesh(in,V,F);

  //   ~~~~~ STEP 1: remesh (copied and pasted from botsch-kobbelt-remesher-libigl/remeshmesh.cpp) ~~~~~

    bool project = false;
    int iterations = 10;
    double h = 0.05;
    
    Eigen::VectorXd target = Eigen::VectorXd::Constant(V.rows(),h);
    Eigen::VectorXi feature;
	feature.resize(0);
    remesh_botsch(V,F,target,iterations, feature, project); // V, F, and target are taken 
  // all_boundary_loop(F);

  //  ~~~~~ STEP 2: VIEW MESH ~~~~~
  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.launch();
}
