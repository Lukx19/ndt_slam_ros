#ifndef GRAPH_SLAM_UK_COVARIANCE_WRAPPER
#define GRAPH_SLAM_UK_COVARIANCE_WRAPPER

#include <graph_slam_uk/utils/eigen_tools.h>
namespace slamuk
{
class CovarianceWrapper
{
public:
  typedef eigt::transform2d_t<double> Transform;
  explicit CovarianceWrapper(const Eigen::Matrix3d &covar) : covar_(covar)
  {
  }
  CovarianceWrapper() : covar_(Eigen::Matrix3d::Identity())
  {
  }

  void addToCovar(const Eigen::Matrix3d &other_covar, const Transform &trans)
  {
    double angle = std::atan2(trans.rotation()(1, 0), trans.rotation()(0, 0));
    double co = std::cos(angle);
    double si = std::sin(angle);
    covar_ =
        computeJacc1(co, si, trans.translation()) * this->covar_ *
            computeJacc1(co, si, trans.translation()).transpose() +
        computeJacc2(co, si) * other_covar * computeJacc2(co, si).transpose();
  }

public:
  Eigen::Matrix3d covar_;

private:
  Eigen::Matrix3d computeJacc1(double co, double si,
                               const Eigen::Vector2d &t_j) const
  {
    Eigen::Matrix3d res;
    res << 1, 0, -si *t_j(0) - co * t_j(1), 0, 1, co * t_j(0) - si * t_j(1), 0,
        0, 1;
    return res;
  }
  Eigen::Matrix3d computeJacc2(double co, double si) const
  {
    Eigen::Matrix3d res;
    res << co, -si, 0, si, co, 0, 0, 0, 1;
    return res;
  }
};
}

#endif
