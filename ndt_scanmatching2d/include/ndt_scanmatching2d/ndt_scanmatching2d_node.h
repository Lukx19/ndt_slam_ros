#ifndef NDT_SCANMATCHING2D_NODE
#define NDT_SCANMATCHING2D_NODE

#include <string>
#include <iostream>
#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <message_filters/synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <laser_geometry/laser_geometry.h>
#include <sensor_msgs/PointCloud2.h>

#include <tf/transform_datatypes.h>
#include <tf/transform_broadcaster.h>

#include <pcl_conversions/pcl_conversions.h>
#include <pcl_ros/transforms.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <Eigen/Dense>
#include <ndt_scanmatching2d/ndt2d.h>
#include <dynamic_slam_utils/eigen_tools.h>

class NdtScanmatching2d
{
  enum ScanmatchingTypes { NDT2D };
  typedef message_filters::sync_policies::ApproximateTime<
      nav_msgs::Odometry, sensor_msgs::LaserScan> ImuSyncPolicy;
  typedef message_filters::Subscriber<nav_msgs::Odometry> odom_sub_t;
  typedef message_filters::Subscriber<sensor_msgs::LaserScan> laser_sub_t;
  typedef Eigen::Vector3d pose_t;
  typedef pcl::PointCloud<pcl::PointXYZ> pcl_t;

public:
  NdtScanmatching2d(ros::NodeHandle &n, ros::NodeHandle &n_private);
  void start();

private:
  ros::NodeHandle nh_;
  ros::NodeHandle nh_private_;
  pcl::NormalDistributionsTransform2DEx<pcl::PointXYZ, pcl::PointXYZ>
      ndt2d_matcher_;

  pcl_t::Ptr old_scan_;
  pose_t old_odom_;
  pose_t old_position_;
  bool is_ready_;
  eigt::transform2d_t<double> tf_trans_;
  // parameters from launch file
  std::string odom_frame_;
  std::string pose_frame_;
  std::string robot_base_frame_;
  std::string tf_prefix_;
  std::string odom_topic_;
  std::string pose_pub_topic_;
  std::string laser_topic_;
  float max_range_;
  size_t resolution_;
  size_t layers_;
  ScanmatchingTypes mode_;
  double min_angle_;
  double min_displacement_;

  laser_geometry::LaserProjection projector_;
  tf::TransformListener tf_list_;
  tf::TransformBroadcaster tf_broadcast_;
  ros::Publisher pose_pub_;
  odom_sub_t odom_sub_;
  laser_sub_t laser_sub_;
  message_filters::Synchronizer<ImuSyncPolicy> msg_sync_;
  uint seq_;
  bool is_initialized;

  void data_cb(const nav_msgs::Odometry::ConstPtr &odom,
               const sensor_msgs::LaserScan::ConstPtr &laser);
  void initParameters();

  bool getTransfNdt2d(eigt::transform2d_t<double> &trans,
                      Eigen::Matrix3d &covar, const pose_t &odom_pose,
                      const pcl_t &points);

  void publishTFTransform(const eigt::transform2d_t<double> &tf_transform,
                          const ros::Time &time);
  void updateTFTrans(const pose_t &odom, const pose_t &calc);


  boost::array<double, 36> convertCovariance(const Eigen::Matrix3d &covar) const;
};

#endif
