#ifndef NDT_GSLAM_COVARIANCE_INVERSE
#define NDT_GSLAM_COVARIANCE_INVERSE
#include <Eigen/Dense>
namespace slamuk
{
/**
 * @brief      Callculate inverse of covariance with singularity check
 *
 * @param[in]  cov           The covariance
 * @param      adjusted_cov  The adjusted covariance
 * @param      icov          The inverse covariance
 *
 * @tparam     Scalar        Scalar of the matrix
 *
 * @return     true if sucesfull
 */
template <typename Scalar>
bool covarInverse(const Eigen::Matrix<Scalar, 3, 3> &cov,
                  Eigen::Matrix<Scalar, 3, 3> &adjusted_cov,
                  Eigen::Matrix<Scalar, 3, 3> &icov)
{
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 3, 3>> solver(cov);
  Eigen::Matrix<Scalar, 3, 3> evecs = solver.eigenvectors().real();
  Eigen::Matrix<Scalar, 3, 3> evals = solver.eigenvalues().asDiagonal();
  double eval_factor = 100;
  if (evals.sum() <= 0) {
    return false;
  } else {
    bool recalc = false;
    // guard against near singular matrices::
    double max_eval = evals.maxCoeff();
    // test if every eigenval is big enough
    for (int i = 0; i < evals.cols(); ++i) {
      if (max_eval > evals(i) * eval_factor) {
        evals(i, i) = max_eval / eval_factor;
        recalc = true;
      }
    }
    if (recalc) {
      adjusted_cov = evecs * evals * (evecs.transpose());
    } else {
      adjusted_cov = cov;
    }
    // compute inverse covariance
    icov = evecs * (evals.inverse()) * (evecs.transpose());
  }
  return true;
}
}

#endif
