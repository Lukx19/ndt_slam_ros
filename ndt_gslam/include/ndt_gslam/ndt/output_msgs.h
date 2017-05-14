#ifndef NDT_GSLAM_NDT_OUTPUT_MSGS
#define NDT_GSLAM_NDT_OUTPUT_MSGS

#include <Eigen/Dense>

namespace slamuk
{
struct OccupancyGrid {
  // origin-----x---->
  // |
  // |
  // y
  // |
  // |
  // v
  //
  /**
  x, y, theta of first node in cell array (left top
  cell in grid) in world's coordinate system
  */
  Eigen::Vector3f origin_;
  size_t width_;
  size_t height_;
  float resolution_;  // [m/cell]
  std::vector<int8_t> cells_;
};

struct NDTCellMsg {
  Eigen::Vector3d mean_;  // x,y,z
  Eigen::Matrix3d cov_;
  double occupancy_;
  long points_;
};

struct NDTGridMsg {
  Eigen::Vector3d origin_;  // x, y, theta of first node in cell array (left top
                            // cell in grid)
  Eigen::Vector3d size_;    // size in x, ,y ,z dimmention
  Eigen::Vector3d cell_sizes_;
  std::vector<NDTCellMsg> cells_;
};
}

#endif
