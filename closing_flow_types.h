#pragma once

struct ClosingFlowParams {
    int maxiter = 20;
    double dt = 0.07;
    double bd = 1.0 / 1.0;

    /// Target edge length for Botsch-Kobbelt remeshing. If <= 0, uses
    double h = 0.85;
    double avg_edge = 0.0;

    /// Remesh iterations per flow step
    int remesh_iterations = 1;

    bool opening = false;
    bool always_recompute = false;
    double tol = 1e-5;
    bool quadric_curvature = false;
    bool verbose = true; 
    Eigen::VectorXi selection;
};