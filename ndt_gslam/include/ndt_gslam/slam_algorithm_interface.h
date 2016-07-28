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
  getOccupancyGrid(const std::string &world_frame_id) const = 0;

  /**
   * @brief      Gets the graph serialized.
   *
   * @param[in]  world_frame_id  The world frame identifier
   *
   * @return     The graph serialized.
   */
  virtual visualization_msgs::MarkerArray
  getGraphSerialized(const std::string &world_frame_id) const = 0;

  /**
   * @brief      Gets the ndt map.
   *
   * @param[in]  world_frame_id  The world frame identifier
   *
   * @return     The ndt map.
   */
  virtual ndt_gslam::NDTMapMsg
  getNDTMap(const std::string &world_frame_id) const = 0;

  /**
   * @brief      Gets the pcl map.
   *
   * @param[in]  world_frame_id  The world frame identifier
   *
   * @return     The pcl map.
   */
  virtual typename PointCloud::Ptr
  getPclMap(const std::string &world_frame_id) const = 0;

  /**
   * @brief      Gets the pcl map 2.
   *
   * @param[in]  world_frame_id  The world frame identifier
   *
   * @return     The pcl map 2.
   */
  virtual typename PointCloud::Ptr
  getPclMap2(const std::string &world_frame_id) const = 0;

  // parameters

  /**
   * @brief      Sets the move window radius.
   *
   * @param[in]  radius  The radius
   */
  virtual void setRunWindowRadius(float radius) = 0;

  /**
   * @brief      Sets the generation distance.
   *
   * @param[in]  distance  The distance
   */
  virtual void setGenerationDistance(float distance) = 0;

  /**
   * @brief      Sets the loop closure maximum distance.
   *
   * @param[in]  dist  The distance
   */
  virtual void setLoopClosureMaxDist(float dist) = 0;

  /**
   * @brief      Sets the loop closure minimum distance.
   *
   * @param[in]  dist  The distance
   */
  virtual void setLoopClosureMinDist(float dist) = 0;

  /**
   * @brief      Sets the loop closure score threshold.
   *
   * @param[in]  score  The score
   */
  virtual void setLoopClosureScoreThreshold(float score) = 0;
};
}
#endif
