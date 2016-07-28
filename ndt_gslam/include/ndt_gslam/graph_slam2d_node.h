#ifndef NDT_GSLAM_GRAPH_SLAM2D_NODE
#define NDT_GSLAM_GRAPH_SLAM2D_NODE

#include <ros/ros.h>

#include <geometry_msgs/Pose.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/LaserScan.h>
#include <sensor_msgs/PointCloud.h>
#include <visualization_msgs/MarkerArray.h>

#include <message_filters/subscriber.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <message_filters/synchronizer.h>
#include <message_filters/time_synchronizer.h>

#include <laser_geometry/laser_geometry.h>

#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_listener.h>
#include <tf_conversions/tf_eigen.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl_ros/point_cloud.h>
#include <pcl_ros/transforms.h>
//#include <pcl/filters/voxel_grid.h>

#include <ndt_gslam/NDTMapMsg.h>
#include <ndt_gslam/slam_algorithm_interface.h>
#include <ndt_gslam/utils/eigen_tools.h>

namespace slamuk
{
class GraphSlamNode
{
  typedef message_filters::sync_policies::ApproximateTime<
      nav_msgs::Odometry, sensor_msgs::LaserScan>
      OdomSyncPolicy;
  typedef message_filters::sync_policies::ApproximateTime<
      geometry_msgs::PoseWithCovarianceStamped, sensor_msgs::LaserScan>
      PoseSyncPolicy;
  typedef message_filters::Subscriber<nav_msgs::Odometry> odom_sub_t;
  typedef message_filters::Subscriber<sensor_msgs::LaserScan> laser_sub_t;
  typedef message_filters::Subscriber<geometry_msgs::PoseWithCovarianceStamped>
      pose_sub_t;
  typedef Eigen::Vector3d pose_t;
  typedef pcl::PointCloud<pcl::PointXYZ> Pcl;
  typedef Pcl::Ptr pcl_ptr_t;
  typedef eigt::transform2d_t<double> transform_t;

public:
  GraphSlamNode(ros::NodeHandle &n, ros::NodeHandle &n_private,
                ISlamAlgorithm &alg);
  void start();

private:
  ros::NodeHandle nh_;
  ros::NodeHandle nh_private_;
  ISlamAlgorithm *algorithm_;
  pose_t last_odom_;
  pose_t last_pose_;

  // ros interface utils
  laser_geometry::LaserProjection projector_;
  tf::TransformListener tf_list_;
  tf::TransformBroadcaster tf_broadcast_;
  odom_sub_t odom_sub_;
  pose_sub_t pose_sub_;
  laser_sub_t laser_sub_;
  ros::Subscriber laser_only_sub_;
  message_filters::Synchronizer<OdomSyncPolicy> odom_sync_;
  message_filters::Synchronizer<PoseSyncPolicy> pose_sync_;
  ros::Publisher occ_map_pub_;
  ros::Publisher win_map_pub_;
  ros::Publisher graph_pub_;
  ros::Publisher win_pcl_pub_;
  ros::Publisher map_pcl_pub_;
  tf::Transform world_tf_trans_;
  uint seq_;
  bool is_ready_;

  // parameters from launch file***********
  std::string fixed_frame_;
  std::string robot_base_frame_;
  std::string odom_frame_;
  std::string tf_prefix_;
  std::string laser_topic_;
  std::string subscribe_mode_;
  // graph serializer info
  bool serialize_graph;
  std::string graph_pub_topic_;

  void initParameters();

  void odom_cb(const nav_msgs::Odometry::ConstPtr &odom,
               const sensor_msgs::LaserScan::ConstPtr &laser);
  void pose_cb(const geometry_msgs::PoseWithCovarianceStamped::ConstPtr &pose,
               const sensor_msgs::LaserScan::ConstPtr &laser);
  void laser_cb(const sensor_msgs::LaserScan::ConstPtr &laser);

  bool prepareAllData(const ros::Time &time_stamp,
                      const sensor_msgs::LaserScan::ConstPtr &laser,
                      pose_t &pose, pcl_ptr_t &pcl);
  bool preparePoseData(const ros::Time &time_stamp,
                       const std::string &laser_frame_id, pose_t &pose,
                       tf::StampedTransform &tf_base);
  bool prepareLaserData(const sensor_msgs::LaserScan::ConstPtr &laser,
                        const tf::StampedTransform &tf_base, pcl_ptr_t &pcl);

  void doAlgorithm(const ros::Time &time_stamp, const pose_t &pose,
                   pcl_ptr_t &pcl, const Eigen::Matrix3d &covar);
  void publishTF(const ros::Time &time);
  void updateTFTrans(const pose_t &odom_pose, const pose_t &calc_pose);
  Eigen::Matrix3d arrayToMatrix(const boost::array<double, 36> &array) const;

  bool movedEnough(const eigt::transform2d_t<double> &trans) const;
  tf::Transform eigenPoseToTF(const Eigen::Vector3d &pose) const;
};
}

#endif
