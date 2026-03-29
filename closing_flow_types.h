#pragma once

struct ClosingFlowParams {
    int maxiter = 1;
    double dt = 0.001;
    double bd = 1.0 / 0.4;

    /// Target edge length for Botsch-Kobbelt remeshing. If <= 0, uses
    /// 0.5 * avg_edge_length (adaptive to mesh scale)
    double h = 0.01;

    /// Remesh iterations per flow step
    int remesh_iterations = 10;

    bool opening = false;
    bool always_recompute = false;
    double tol = 1e-5;
};
