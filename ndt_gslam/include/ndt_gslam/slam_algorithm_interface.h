#ifndef NDT_GSLAM_SLAM_ALGORITHM_INTERFACE
#define NDT_GSLAM_SLAM_ALGORITHM_INTERFACE

#include <nav_msgs/OccupancyGrid.h>
#include <ndt_gslam/NDTMapMsg.h>
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

  /**
   * @brief      This function updates state of the whole algorithm.
   *
   * @param[in]  motion       The motion from odom if available
   * @param[in]  covariance   The covariance
   * @param[in]  pcl          The current scan
   * @param[in]  update_time  The update time
   *
   * @return     current position in the world frame
   */
  virtual Pose update(const Transform &motion, const Covar &covariance,
                      const PointCloud &pcl, const ros::Time &update_time) = 0;

  /**
   * @brief      Gets the occupancy grid.
   *
   * @param[in]  world_frame_id  The world frame identifier
   *
   * @return     The occupancy grid.
   */
  virtual nav_msgs::OccupancyGrid
  getOccupancyGrid(const std::string &world_frame_id) = 0;

  /**
   * @brief      Gets the graph serialized.
   *
   * @param[in]  world_frame_id  The world frame identifier
   *
   * @return     The graph serialized.
   */
  virtual visualization_msgs::MarkerArray
  getGraphSerialized(const std::string &world_frame_id) = 0;

  /**
   * @brief      Gets the ndt map.
   *
   * @param[in]  world_frame_id  The world frame identifier
   *
   * @return     The ndt map.
   */
  virtual ndt_gslam::NDTMapMsg getNDTMap(const std::string &world_frame_id) = 0;
};
}
#endif
