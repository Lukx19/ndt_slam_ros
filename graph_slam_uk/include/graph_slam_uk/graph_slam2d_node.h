#ifndef GRAPH_SLAM_UK_GRAPH_SLAM2D_NODE
#define GRAPH_SLAM_UK_GRAPH_SLAM2D_NODE

#include <ros/ros.h>

#include <sensor_msgs/LaserScan.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/PointCloud.h>
#include <visualization_msgs/MarkerArray.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>

#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <message_filters/synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>

#include <laser_geometry/laser_geometry.h>

#include <tf/transform_datatypes.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_listener.h>
#include <tf_conversions/tf_eigen.h>

#include <pcl_conversions/pcl_conversions.h>
#include <pcl_ros/transforms.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
//#include <pcl/filters/voxel_grid.h>

#include <graph_slam_uk/slam_optimizer/graph_slam_interfaces.h>
#include <graph_slam_uk/utils/eigen_tools.h>

namespace slamuk
{
template <typename T>
class GraphSlamNode
{
  typedef message_filters::sync_policies::ApproximateTime<
      nav_msgs::Odometry, sensor_msgs::LaserScan> OdomSyncPolicy;
  typedef message_filters::sync_policies::ApproximateTime<
      geometry_msgs::PoseWithCovarianceStamped, sensor_msgs::LaserScan>
      PoseSyncPolicy;
  typedef message_filters::Subscriber<nav_msgs::Odometry> odom_sub_t;
  typedef message_filters::Subscriber<sensor_msgs::LaserScan> laser_sub_t;
  typedef message_filters::Subscriber<geometry_msgs::PoseWithCovarianceStamped>
      pose_sub_t;
  typedef Eigen::Vector3d pose_t;
  typedef pcl::PointCloud<pcl::PointXYZ> pcl_t;
  typedef pcl_t::Ptr pcl_ptr_t;
  typedef eigt::transform2d_t<double> transform_t;

public:
  GraphSlamNode(ros::NodeHandle &n, ros::NodeHandle &n_private,
                IGraphOptimalizer2d<T> &opt_engine, IScanmatcher2d &scan_match);
  void start();

private:
  ros::NodeHandle nh_;
  ros::NodeHandle nh_private_;
  IGraphOptimalizer2d<T> *opt_engine_;
  IScanmatcher2d *scan_match_;
  std::vector<size_t> scans_;
  pose_t last_odom_;
  pcl_ptr_t last_scan_;
  pcl_ptr_t all_scans_;

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
  ros::Publisher pcl_map_pub_;
  ros::Publisher graph_pub_;
  transform_t world_tf_trans_;
  uint seq_;
  bool is_ready_;

  // parameters from launch file***********
  std::string fixed_frame_;
  std::string robot_base_frame_;
  std::string odom_frame_;
  std::string tf_prefix_;
  std::string odom_topic_;
  std::string pose_topic_;
  std::string pclmap_pub_topic_;
  std::string laser_topic_;
  std::string subscribe_mode_;
  // graph serializer info
  bool serialize_graph;
  std::string graph_pub_topic_;
  // // scanmatching atributes
  float max_range_;
  double grid_step_;
  double min_rotation_;
  double min_displacement_;
  // // optimalizer atributes
  size_t iterations_;
  double epsilon_err_;

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
  void updateTFTrans(pose_t odom_pose, pose_t calc_pose);
  Eigen::Matrix3d arrayToMatrix(const boost::array<double, 36> &array) const;

  bool movedEnough(const eigt::transform2d_t<double> &trans) const;
  void saveDotGraph();
};
// ***************************IMPLEMENTATION *************************
template <typename T>
GraphSlamNode<T>::GraphSlamNode(ros::NodeHandle &n, ros::NodeHandle &n_private,
                                IGraphOptimalizer2d<T> &opt_engine,
                                IScanmatcher2d &scan_match)
  : nh_(n)
  , nh_private_(n_private)
  , opt_engine_(&opt_engine)
  , scan_match_(&scan_match)
  , all_scans_(new pcl_t())
  , odom_sync_(OdomSyncPolicy(10))
  , pose_sync_(PoseSyncPolicy(10))
  , world_tf_trans_(transform_t::Identity())
  , seq_(0)
  , is_ready_(false)
{
  initParameters();
  ROS_INFO("Graph_slam2d: Launch params initialized.");
  scan_match_->setGridStep(grid_step_);
  scan_match_->setMaxRange(max_range_);
  scan_match_->setTransformationEpsilon(epsilon_err_);

  opt_engine_->setMaxIterations(iterations_);
  opt_engine_->setEuclideanMaxError(epsilon_err_);

  pcl_map_pub_ =
      nh_.advertise<sensor_msgs::PointCloud2>(pclmap_pub_topic_, 5, false);
  graph_pub_ = nh_.advertise<visualization_msgs::MarkerArray>(graph_pub_topic_,
                                                              5, false);
  const int message_cache = 10;
  if (subscribe_mode_ == "ODOM") {
    laser_sub_.subscribe(nh_, laser_topic_, message_cache);
    odom_sub_.subscribe(nh_, odom_topic_, message_cache);
    // sync messages using approximate alghorithm
    odom_sync_.connectInput(odom_sub_, laser_sub_);
    odom_sync_.registerCallback(
        boost::bind(&GraphSlamNode::odom_cb, this, _1, _2));
  } else if (subscribe_mode_ == "POSE") {
    laser_sub_.subscribe(nh_, laser_topic_, message_cache);
    pose_sub_.subscribe(nh_, pose_topic_, message_cache);
    // sync messages using approximate alghorithm
    pose_sync_.connectInput(pose_sub_, laser_sub_);
    pose_sync_.registerCallback(
        boost::bind(&GraphSlamNode::pose_cb, this, _1, _2));
  } else {
    // NON
    laser_only_sub_ = nh_.subscribe<sensor_msgs::LaserScan>(
        laser_topic_, message_cache,
        boost::bind(&GraphSlamNode::laser_cb, this, _1));
  }

  ROS_INFO("Graph_slam2d: Node is initialized.");
}

template <typename T>
void GraphSlamNode<T>::start()
{
  ros::spin();
}

template <typename T>
void GraphSlamNode<T>::initParameters()
{
  // find tf_prefix if exists add it to tf_prefix_ variable
  tf_prefix_ = "";
  std::string tf_prefix_path;
  if (nh_private_.searchParam("tf_prefix", tf_prefix_path)) {
    nh_private_.getParam(tf_prefix_path, tf_prefix_);
  }

  robot_base_frame_ =
      nh_private_.param<std::string>("robot_base_frame_id", "base_link");

  odom_frame_ = nh_private_.param<std::string>("odom_farme_id", "odom");

  fixed_frame_ = nh_private_.param<std::string>("fixed_farme_id", "world");

  odom_topic_ = nh_private_.param<std::string>("odom_topic", "/odom");

  pose_topic_ = nh_private_.param<std::string>("pose_topic", "/odom");

  subscribe_mode_ = nh_private_.param<std::string>("subscribe_mode", "NON");

  pclmap_pub_topic_ =
      nh_private_.param<std::string>("pcl_map_topic", "/pose_calc");

  laser_topic_ = nh_private_.param<std::string>("laser_topic", "/laser");

  grid_step_ = nh_private_.param<double>("grid_step", 0.5);

  max_range_ = static_cast<float>(
      nh_private_.param<double>("maximal_laser_range", 30.0));

  min_rotation_ = nh_private_.param<double>("min_rotated_angle", 0);

  min_displacement_ = nh_private_.param<double>("min_traveled_distance", 0);

  iterations_ =
      static_cast<size_t>(nh_private_.param<int>("optimalizer_iterations", 5));

  epsilon_err_ = nh_private_.param<double>("min_epsilon_change", 0.001);

  serialize_graph = nh_private_.param<bool>("serialize_graph", true);
  graph_pub_topic_ = nh_private_.param<std::string>("graph_pub_topic", "graph");

  if (tf_prefix_ != "") {
    robot_base_frame_ = tf_prefix_ + "/" + robot_base_frame_;
    odom_frame_ = tf_prefix_ + "/" + odom_frame_;
    fixed_frame_ = tf_prefix_ + "/" + fixed_frame_;
  }
}
template <typename T>
void GraphSlamNode<T>::odom_cb(const nav_msgs::Odometry::ConstPtr &odom,
                               const sensor_msgs::LaserScan::ConstPtr &laser)
{
  pose_t new_pose;  // current robot position on the odom_frame
  pcl_ptr_t pcl;    // pcl holding laser scan transformed in base_frame of robot
  ros::Time t = odom->header.stamp;
  if (!prepareAllData(t, laser, new_pose, pcl)) {
    publishTF(t);
    return;
  }
  doAlgorithm(t, new_pose, pcl, arrayToMatrix(odom->pose.covariance));
}

template <typename T>
void GraphSlamNode<T>::pose_cb(
    const geometry_msgs::PoseWithCovarianceStamped::ConstPtr &pose,
    const sensor_msgs::LaserScan::ConstPtr &laser)
{
  pose_t new_pose;  // current robot position on the odom_frame
  pcl_ptr_t pcl;    // pcl holding laser scan transformed in base_frame of robot
  ros::Time t = pose->header.stamp;
  if (!prepareAllData(t, laser, new_pose, pcl)) {
    publishTF(t);
    return;
  }
  doAlgorithm(t, new_pose, pcl, arrayToMatrix(pose->pose.covariance));
}
template <typename T>
void GraphSlamNode<T>::laser_cb(const sensor_msgs::LaserScan::ConstPtr &laser)
{
  pose_t new_pose;  // current robot position on the odom_frame
  pcl_ptr_t pcl;    // pcl holding laser scan transformed in base_frame of robot
  ros::Time t = laser->header.stamp;
  if (!prepareAllData(t, laser, new_pose, pcl)) {
    publishTF(t);
    return;
  }
  Eigen::Matrix3d covar;
  // add here more sophisticated aproximation of covariance in future
  covar << 1e-4, 0, 0, 0, 1e-4, 0, 0, 0, 1e-4;
  doAlgorithm(t, new_pose, pcl, covar);
}

template <typename T>
bool
GraphSlamNode<T>::prepareAllData(const ros::Time &time_stamp,
                                 const sensor_msgs::LaserScan::ConstPtr &laser,
                                 pose_t &pose, pcl_ptr_t &pcl)
{
  tf::StampedTransform tf_base;  // transformation from laser_frame ->
                                 // base_frame
  if (!preparePoseData(time_stamp, laser->header.frame_id, pose, tf_base)) {
    return false;
  }

  if (!prepareLaserData(laser, tf_base, pcl)) {
    return false;
  }
  ROS_INFO("GraphSlam: Messages transformed.");
  return true;
}

template <typename T>
bool GraphSlamNode<T>::preparePoseData(const ros::Time &time_stamp,
                                       const std::string &laser_frame_id,
                                       pose_t &pose,
                                       tf::StampedTransform &tf_base)
{
  // get transformations from TF
  tf::StampedTransform tf_odom;
  try {
    tf_list_.waitForTransform(robot_base_frame_, laser_frame_id, time_stamp,
                              ros::Duration(10.0));
    tf_list_.lookupTransform(robot_base_frame_, laser_frame_id, time_stamp,
                             tf_base);

    tf_list_.waitForTransform(odom_frame_, robot_base_frame_, time_stamp,
                              ros::Duration(10.0));
    tf_list_.lookupTransform(odom_frame_, robot_base_frame_, time_stamp,
                             tf_odom);
  } catch (...) {
    ROS_ERROR_STREAM(
        "GraphSlam: Not sucessfull in retrieving tf tranform laser "
        "-> base -> odom ");
    return false;
  }
  Eigen::Affine3d trans_robot;
  tf::transformTFToEigen(tf_odom, trans_robot);
  pose = eigt::getPoseFromTransform(
      eigt::convertToTransform(trans_robot.matrix()));

  // check if robot moved enough
  auto odom_trans_mat = eigt::transBtwPoses(last_odom_, pose);
  if (!movedEnough(odom_trans_mat)) {
    ROS_INFO("Graph_slam2d:Robot not moved enough-> droping msgs");
    return false;
  }
  return true;
}

template <typename T>
bool GraphSlamNode<T>::prepareLaserData(
    const sensor_msgs::LaserScan::ConstPtr &laser,
    const tf::StampedTransform &tf_base, pcl_ptr_t &pcl)
{
  // project laser message to point cloud class in laser frame_id
  sensor_msgs::PointCloud2 laser_pcl_msg;
  pcl_t laser_pcl;
  pcl_t laser_pcl_base;
  /////////////
  projector_.projectLaser(*laser, laser_pcl_msg);
  pcl::moveFromROSMsg(laser_pcl_msg, laser_pcl);
  // transform point cloud from laser frame_id -> robot base frame
  try {
    pcl_ros::transformPointCloud(laser_pcl, laser_pcl_base, tf_base);
  } catch (tf::TransformException &e) {
    ROS_ERROR_STREAM("GraphSlam: Point cloud transformation not successful"
                     << e.what());
    return false;
  }
  if (laser_pcl_base.size() < 50) {
    ROS_ERROR_STREAM("GraphSlam: Not enough points in laser scan "
                     << laser_pcl_base.size());
    return false;
  }
  // make copy of pointcloud wrapped in shared ptr
  pcl = laser_pcl_base.makeShared();
  return true;
}

template <typename T>
void GraphSlamNode<T>::doAlgorithm(const ros::Time &time_stamp,
                                   const pose_t &pose, pcl_ptr_t &pcl,
                                   const Eigen::Matrix3d &covar)
{
  // prepare initial data in first received msg pair
  if (!is_ready_) {
    (*all_scans_) += (*pcl);
    last_odom_ = pose;
    last_scan_ = pcl;
    scans_.push_back(opt_engine_->addPose(pose, pcl));
    is_ready_ = true;
    publishTF(time_stamp);
    return;
  }
  // adding odometry constrain to optimalizer
  pose_t odom_trans =
      eigt::getPoseFromTransform(eigt::transBtwPoses(last_odom_, pose));
  // auto res_match =
  // scan_match_->match(pcl,last_scan_,eigt::transBtwPoses(last_odom_,
  // pose).matrix());
  // opt_engine_->addLastConstrain(odom_trans, covar.inverse());
  // if(res_match.success_){
  //   ROS_DEBUG("using scanmatched edge");
  //   // adding node to optimalizer;
  //   auto matcher_pose = eigt::transformPose(last_odom_,res_match.transform_);
  //   scans_.push_back(opt_engine_->addPose(matcher_pose, pcl));
  //   opt_engine_->addLastConstrain(eigt::getPoseFromTransform(res_match.transform_),
  //   res_match.inform_);
  //   updateTFTrans(pose,matcher_pose);
  //   last_odom_ = matcher_pose;
  // }else{
  // adding node to optimalizer;
  scans_.push_back(opt_engine_->addPose(pose, pcl));
  opt_engine_->addLastConstrain(odom_trans, covar.inverse());
  last_odom_ = pose;
  // }
  // try to close loop closures

  bool res = false;
  if (opt_engine_->tryLoopClose()) {
    res = opt_engine_->optimalize();
  }
  ROS_INFO("Graph_slam2d:Calculation Finished");

  // update old point cloud or recreate new one
  pcl_t pcl_transformed;
  if (res) {
    all_scans_->clear();
    for (auto id : scans_) {
      auto trans4x4 = eigt::convertFromTransform(
          eigt::getTransFromPose(opt_engine_->getPoseLocation(id)));
      pcl::transformPointCloud(*opt_engine_->getPoseData(id), pcl_transformed,
                               trans4x4);
      (*all_scans_) += pcl_transformed;
    }
  } else {
    size_t last_id = scans_.back();
    auto last_pose = opt_engine_->getPoseLocation(last_id);
    auto trans4x4 =
        eigt::convertFromTransform(eigt::getTransFromPose(last_pose));
    pcl::transformPointCloud((*pcl), pcl_transformed, trans4x4);
    (*all_scans_) += pcl_transformed;
  }

  // PUBLISHING MSGS

  // downsample pointcloud
  // not working currently with c++11
  // pcl::VoxelGrid<pcl::PointXYZ> sor;
  // sor.setInputCloud (all_scans_);
  // sor.setLeafSize (0.5f, 0.5f, 0.5f);
  // pcl_ptr_t filtered_cloud;
  // sor.filter (*filtered_cloud);
  // publish msgs with map
  sensor_msgs::PointCloud2 msg;
  pcl::toROSMsg(*all_scans_, msg);
  msg.header.frame_id = fixed_frame_;
  msg.header.stamp = time_stamp;
  pcl_map_pub_.publish(msg);
  ROS_INFO("Graph_slam2d:PCL Map published.");

  // updateTFTrans(pose,opt_engine_->getPoseLocation(scans_.back()));
  publishTF(time_stamp);

  // publish serialized graph
  if (serialize_graph) {
    graph_pub_.publish(opt_engine_->getGraphSerialized(fixed_frame_));
    // saveDotGraph();
    ROS_INFO("Graph_slam2d:GRAPH Markers published.");
  }
  // prepare data for next iteration
  // last_odom_ = pose;
  last_scan_ = pcl;
  ++seq_;
}

template <typename T>
void GraphSlamNode<T>::publishTF(const ros::Time &time)
{
  // publish map tf transform
  auto w_trans = eigt::getPoseFromTransform(world_tf_trans_);
  tf::Transform transform;
  transform.setOrigin(tf::Vector3(w_trans(0), w_trans(1), 0.0));
  tf::Quaternion q;
  q.setRPY(0, 0, w_trans(2));
  transform.setRotation(q);
  tf_broadcast_.sendTransform(
      tf::StampedTransform(transform, time, fixed_frame_, odom_frame_));
}

template <typename T>
void GraphSlamNode<T>::updateTFTrans(pose_t odom, pose_t slam)
{
  auto odom_w_frame = eigt::transformPose(odom, world_tf_trans_);
  world_tf_trans_ = eigt::transBtwPoses(odom_w_frame, slam) * world_tf_trans_;
}

template <typename T>
Eigen::Matrix3d
GraphSlamNode<T>::arrayToMatrix(const boost::array<double, 36> &array) const
{
  Eigen::Matrix3d covar;
  covar << array[0], array[1], array[5], array[6], array[7], array[11],
      array[30], array[31], array[35];
  return covar;
}

template <typename T>
bool
GraphSlamNode<T>::movedEnough(const eigt::transform2d_t<double> &trans) const
{
  if (eigt::getAngle(trans) < min_rotation_ &&
      eigt::getDisplacement(trans) < min_displacement_)
    return false;
  return true;
}

template <typename T>
void GraphSlamNode<T>::saveDotGraph()
{
  std::ofstream out;
  out.open("graph" + std::to_string(ros::Time::now().toSec()) + ".dot");
  opt_engine_->getGraphSerialized(out);
  out.close();
}
}

#endif
