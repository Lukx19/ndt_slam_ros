#ifndef NDT_GRID2D_NDT_GRID2D_INTERFACE
#define NDT_GRID2D_NDT_GRID2D_INTERFACE
#include <Eigen/Dense>
#include <pcl/point_cloud.h>

struct OccupancyGrid
{
  Eigen::Vector3d origin_;
  size_t width_;
  size_t height_;
  float resolution_;  // [m/cell]
  std::vector<int8_t> cells_;
};

template <typename PointType>
class INDTGrid2D
{
public:
  typedef Eigen::Vector3d Pose;
  typedef pcl::PointCloud<PointType> PointCloud;
  typedef INDTGrid2D<PointType> SelfType;

  void initialize(const PointCloud &pcl) = 0;

  void mergeIn(const PointCloud &pcl, const Pose &origin,
               bool resize = true) = 0;
  void mergeIn(const SelfType &grid, bool transform = true,
               bool resize = true) = 0;
  void mergeInTraced(const PointCloud &pcl, const Pose &origin,
                     bool resize = true) = 0;
  void mergeInTraced(const SelfType &grid, bool transform,
                     bool resize = true) = 0;

  void enlarge(float left, float down, float right, float up) = 0;
  void move(const Eigen::Vector2d &translation) = 0;
  void transform(const Eigen::Matrix3d &transform) = 0;
  OccupancyGrid createOccupancyGrid() const = 0;
};
#endif
