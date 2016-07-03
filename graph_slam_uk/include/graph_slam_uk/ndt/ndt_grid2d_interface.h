#ifndef GRAPH_SLAM_UK_NDT_GRID2D_INTERFACE
#define GRAPH_SLAM_UK_NDT_GRID2D_INTERFACE

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

class INDTMappable
{
public:
  virtual ~INDTMappable()
  {
  }
  typedef Eigen::Vector3d Pose;
  typedef INDTMappable SelfType;
  virtual void mergeIn(const SelfType &grid, bool transform = true,
                       bool resize = true) = 0;
  virtual void mergeInTraced(const SelfType &grid, bool transform,
                             bool resize = true) = 0;

  virtual double getTimestamp() const = 0;

  virtual void setTimestamp(double timestamp) = 0;

  virtual bool operator<(const SelfType &other) const = 0;
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

  void setCellSize(float cell_size) = 0;

  Eigen::Vector3d getOrigin() const = 0;

  void setOrigin(const Eigen::Vector3d &origin) = 0;

  double getTimestamp() const = 0;

  void setTimestamp(double timestamp) = 0;

  bool operator<(const INDTGrid2D &other) const = 0;
};
#endif
