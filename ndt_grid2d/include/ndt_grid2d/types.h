#ifndef NDT_GRID2D_TYPES
#define NDT_GRID2D_TYPES

#include <vector>
#include <Eigen/Dense>

namespace slamuk
{
struct OccupancyGrid
{
  Eigen::Vector3d origin_;
  size_t width_;
  size_t height_;
  float resolution_;  // [m/cell]
  std::vector<int8_t> cells_;
};
} // end of namespace slamuk

#endif

