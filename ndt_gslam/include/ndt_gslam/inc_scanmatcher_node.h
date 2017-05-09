#ifndef NDT_GSLAM_INC_SCANMATCHER_NODE
#define NDT_GSLAM_INC_SCANMATCHER_NODE

#include <laser_geometry/laser_geometry.h>
#include <nav_msgs/Odometry.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl_ros/point_cloud.h>
#include <pcl_ros/transforms.h>
#include <ros/console.h>
#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <sensor_msgs/PointCloud.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_listener.h>
#include <tf_conversions/tf_eigen.h>
#include <Eigen/Dense>

#include <pcl/common/transforms.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <ndt_gslam/ndt/ndt_cell.h>
#include <ndt_gslam/ndt/ndt_grid2d.h>
#include <ndt_gslam/registration/d2d_ndt2d.h>
#include <ndt_gslam/registration/ndt2d.h>
#include <ndt_gslam/utils/eigen_tools.h>
#include <ndt_gslam/utils/msgs_conversions.h>
#include <ndt_gslam/utils/tf_conversions.h>

class IncScanmatcherNode
{
public:
  typedef Eigen::Vector3d Pose;
  typedef pcl::PointCloud<pcl::PointXYZ> Pcl;
  typedef eigt::transform2d_t<double> Transform;
  typedef pcl::PointXYZ PointType;
  typedef slamuk::NDTCell CellType;
  typedef slamuk::NDTGrid2D<CellType, PointType> FrameType;
  typedef FrameType::Ptr FrameTypePtr;

public:
  IncScanmatcherNode(ros::NodeHandle &n, ros::NodeHandle &n_private);
  void start();

private:
  ros::NodeHandle nh_;
  ros::NodeHandle nh_private_;
  Transform last_odom_;
  Transform position_;

  Transform unused_trans_;
  Transform unused_odom_;
  Transform position_cumul_;
  // ros interface utils
  laser_geometry::LaserProjection projector_;
  tf::TransformListener tf_list_;
  tf::TransformBroadcaster tf_broadcast_;
  ros::Subscriber laser_sub_;
  ros::Publisher odom_pub_;
  ros::Publisher pcl_pub_;
  uint seq_;

  // parameters from launch file***********
  std::string fixed_frame_;
  std::string robot_base_frame_;
  std::string odom_frame_;
  std::string tf_prefix_;
  std::string laser_topic_;
  float win_radius_;
  float cell_size_;
  double min_rotation_;
  double min_distance_;

  bool initialized_;
  FrameTypePtr running_window_;
  pcl::NormalDistributionsTransform2DEx<PointType, PointType, CellType>
      inc_matcher_;
  //  pcl::D2DNormalDistributionsTransform2D<PointType, PointType, CellType>
  //      inc_matcher_;

  void initParameters();

  void laserCb(const sensor_msgs::LaserScan::ConstPtr &laser);
  Pcl::Ptr prepareLaserData(const sensor_msgs::LaserScan::ConstPtr &laser,
                            const tf::StampedTransform &tf_base);
  bool preparePoseData(const ros::Time &time_stamp,
                       const std::string &laser_frame_id,
                       Transform &odom_position, tf::StampedTransform &tf_base);
  Transform calcScanMovement(const Pcl::Ptr &pcl);

  void publishOdometry(const ros::Time &time);
  void publishPcl(const ros::Time &time);
  bool movedEnough(const Transform &trans) const;
};

#endif  // NDT_GSLAM_INC_SCANMATCHER_NODE
