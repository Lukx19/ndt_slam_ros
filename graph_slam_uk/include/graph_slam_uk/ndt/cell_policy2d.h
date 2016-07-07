#ifndef GRAPH_SLAM_UK_CELL_POLICY2D
#define GRAPH_SLAM_UK_CELL_POLICY2D

#include <Eigen/Dense>

namespace slamuk
{
class CellPolicy2d
{
public:
  typedef Eigen::Vector2d Vector;
  typedef Eigen::Matrix2d Matrix;
  typedef Eigen::Transform<double, 2, Eigen::TransformTraits::Affine> Transform;

  static const size_t dim_ = 2;
  static const size_t max_points_ = 1e9;
  static constexpr float max_occupancy_ = 255.0f;
  static constexpr float sensor_noise_ = 0.01f;
  // std::log(0.6f / (1.0f - 0.6f));
  static constexpr double log_like_occ_ = 0.405465108;
};
}

#endif