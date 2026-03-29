#include "closing_flow.h"

#include <algorithm>
#include <chrono>
#include <iostream>

#include <Eigen/Sparse>

#include <igl/adjacency_matrix.h>
#include <igl/avg_edge_length.h>
#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>

namespace cf = closing_flow_detail;

bool closing_flow(
    const Eigen::MatrixXd &V_in,
    const Eigen::MatrixXi &F_in,
    const ClosingFlowParams &params,
    Eigen::MatrixXd &V_out,
    Eigen::MatrixXi &F_out)
{
    if (F_in.cols() != 3 || V_in.cols() != 3)
        return false;

    const auto t_flow_start = std::chrono::steady_clock::now();
    double remesh_seconds_total = 0;

    Eigen::MatrixXd Vfull = V_in;
    Eigen::MatrixXi Ffull = F_in;

    std::cerr << "Shape of Vfull: (" << Vfull.rows() << ", " << Vfull.cols()
              << ")\n";
    std::cerr << "Shape of Ffull: (" << Ffull.rows() << ", " << Ffull.cols()
              << ")\n";

    bool recompute = true;

    // !!!!!!!!!!!!! NOTE: MAXITER IS HARDCODED TO 1 !!!!!!!!!!!!!
    // When we actually run the algorithm it will be like idk 300 or something big
    for (int iter=0; iter<params.maxiter; iter++) {
        
        Eigen::MatrixXd Vprev = Vfull; // (n, 3)
        Eigen::MatrixXi Fprev = Ffull; // (m, 3)

        const int nVerts = Vfull.rows();
        Eigen::VectorXi moving(nVerts);
        Eigen::VectorXi frozen(nVerts);
        
        if (iter == 0 || recompute) {
            // Compute A, sparse matrix with adjacency matrix plus identity for self-loops
            Eigen::SparseMatrix<int> A;
            igl::adjacency_matrix(Ffull, A); // (n, n)
            for (int i = 0; i < Vfull.rows(); i++) {
                A.coeffRef(i, i) = 1; // (n, n)
            }

            // mass matrix M (use voronoi mass matrix - unspecified in paper?)
            Eigen::SparseMatrix<double> M_sparse;
            igl::massmatrix(Vfull, Ffull, igl::MASSMATRIX_TYPE_VORONOI, M_sparse); // (n, n)
            
            // diagonal elements and clip, vectorize M_spares to get M
            Eigen::VectorXd M(Vfull.rows()); // (n, 1)
            for (int i = 0; i < Vfull.rows(); i++) {
                M(i) = std::max(M_sparse.coeff(i, i), 1e-8);
            }

            // Gaussian curvature K
            Eigen::VectorXd K; // (n, 1)
            igl::gaussian_curvature(Vfull, Ffull, K);

            // Discrete mean curvature H
            Eigen::VectorXd H = cf::discrete_mean_curvature(Vfull, Ffull); // (n, 1)

            // Principal curvatures
            Eigen::VectorXd squaredTerm = (H.array().square() - M.array() * K.array()).max(0.0);
            Eigen::VectorXd sqrtDisc = squaredTerm.array().sqrt();
            Eigen::VectorXd k_max = H + sqrtDisc;
            Eigen::VectorXd k_min = H - sqrtDisc;

            // Select curvature based on opening vs closing
            Eigen::VectorXd curv; // (n, 1) the curvature of each vertex
            if (params.opening) { // opening: NEGATIVE of maximum curvature (so single comparison later is correct for both opening and closing)
                curv = (-k_max.array() / M.array()).matrix(); 
            } else { // closing: minimum curvature
                curv = (k_min.array() / M.array()).matrix();
            }

            // Active vertices where curv < -bd
            // Closing: k_min < -bd
            // Opening: k_max > bd <-> -k_max < -bd
            // allocate at most #V entries then shrink at end so we 
            // don't have to loop through all vertices to form frozen 
            // set dtype Eigen::VectorXi for remeshing
            Eigen::Index nmov = 0, nfroz = 0;
            for (Eigen::Index i = 0; i < nVerts; ++i) {
                if (curv(i) < -params.bd) {
                    moving(nmov++) = static_cast<int>(i);
                } else {
                    frozen(nfroz++) = static_cast<int>(i);
                }
            }
            moving.conservativeResize(nmov);
            frozen.conservativeResize(nfroz);
            std::cerr << "num active vertices: " << moving.size() << "\n";
            std::cerr << "num frozen vertices: " << frozen.size() << "\n";
            std::cerr << "active vertices + frozen vertices: " << moving.size() + frozen.size()
                      << "\n";
        }



        
        // Step 4: Local remeshing (1 iteration)
        if (params.remesh_iterations > 0) {
            double h_target = params.h;
            if (params.h <= 0) {
                double avg = igl::avg_edge_length(Vfull, Ffull);
                
                std::cerr << "avg edge length (avg): " << avg << "\n";

                h_target = 0.1 * avg;
                // h_target = avg;
            }
            std::cerr << "target edge length (h_target): " << h_target << "\n";

            // empty vector for fixed vertices (no frozen vertices)
            // Eigen::VectorXi fixed_vertices;
            // fixed_vertices.resize(0); 

            // expects the indices of the frozen vertices, used as a row index into Vfull
            Eigen::VectorXi fixed_vertices(4); 
            fixed_vertices << 0, 1, 2, 3;

            // Construct target edge length vector
            const int numV = static_cast<int>(Vfull.rows());
            Eigen::VectorXd targetVec = Eigen::VectorXd::Constant(numV, h_target);
            std::cerr << "Starting remeshing\n";
            const auto t_remesh_start = std::chrono::steady_clock::now();
            remesh_botsch(Vfull, Ffull, targetVec, params.remesh_iterations, frozen, true);
            remesh_seconds_total +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                              t_remesh_start)
                    .count();

            // can test with cube and use fixed_vertices
            // keeps one face of the cube fixed
            // remesh_botsch(Vfull, Ffull, targetVec, params.remesh_iterations, fixed_vertices, true);

            std::cerr << "Finished remeshing\n";
        }
    }

    const double flow_seconds_total =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t_flow_start)
            .count();
    std::cerr << "closing_flow total time: " << flow_seconds_total << " s\n";
    std::cerr << "remeshing total time: " << remesh_seconds_total << " s\n";

    V_out = std::move(Vfull);
    F_out = std::move(Ffull);
    return true;
}
