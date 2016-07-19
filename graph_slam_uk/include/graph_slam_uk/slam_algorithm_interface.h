#ifndef GRAPH_SLAM_UK_SLAM_ALGORITHM_INTERFACE
#define GRAPH_SLAM_UK_SLAM_ALGORITHM_INTERFACE

#include <graph_slam_uk/NDTMapMsg.h>
#include <nav_msgs/OccupancyGrid.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <ros/ros.h>
#include <visualization_msgs/MarkerArray.h>
#include <Eigen/Dense>

namespace slamuk
{
class ISlamAlgorithm
{
public:
  typedef Eigen::Transform<double, 2, Eigen::TransformTraits::Affine> Transform;
  typedef Eigen::Matrix3d Covar;
  typedef pcl::PointXYZ PointType;
  typedef pcl::PointCloud<PointType> PointCloud;
  typedef Eigen::Vector3d Pose;

  virtual ~ISlamAlgorithm()
  {
  }

  virtual Pose update(const Transform &motion, const Covar &covariance,
                      const PointCloud &pcl, const ros::Time &update_time) = 0;
  virtual nav_msgs::OccupancyGrid
  getOccupancyGrid(std::string &world_frame_id) const = 0;

  virtual visualization_msgs::MarkerArray
  getGraphSerialized(std::string world_frame_id) const = 0;

  virtual graph_slam_uk::NDTMapMsg
  getNDTMap(std::string &world_frame_id) const = 0;

  virtual typename PointCloud::Ptr
  getPclMap(std::string &world_frame_id) const = 0;

  virtual typename PointCloud::Ptr
  getPclMap2(std::string &world_frame_id) const = 0;

  // parameters
  virtual void setRunWindowRadius(float radius) = 0;
};
}
#endif
