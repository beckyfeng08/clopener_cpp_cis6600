#pragma once

#include <Eigen/Sparse>

namespace closing_flow_detail {

Eigen::MatrixXd dihedral_angles(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    Eigen::MatrixXi& TT_out);

Eigen::VectorXd discrete_mean_curvature(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F);

} // namespace closing_flow_detail
