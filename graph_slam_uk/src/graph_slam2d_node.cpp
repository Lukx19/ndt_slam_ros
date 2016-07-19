#include <graph_slam_uk/graph_slam2d_node.h>
#include <graph_slam_uk/ndt_slam_algorithm.h>

using namespace slamuk;

GraphSlamNode::GraphSlamNode(ros::NodeHandle &n, ros::NodeHandle &n_private,
                             ISlamAlgorithm &alg)
  : nh_(n)
  , nh_private_(n_private)
  , algorithm_(&alg)
  , last_odom_(0, 0, 0)
  , last_pose_(pose_t::Zero())
  , odom_sync_(OdomSyncPolicy(10))
  , pose_sync_(PoseSyncPolicy(10))
  // , world_tf_trans_()
  , seq_(0)
  , is_ready_(false)
{
  initParameters();
  ROS_INFO("Graph_slam2d: Launch params initialized.");

  occ_map_pub_ =
      nh_.advertise<nav_msgs::OccupancyGrid>(map_pub_topic_, 5, false);
  graph_pub_ = nh_.advertise<visualization_msgs::MarkerArray>(graph_pub_topic_,
                                                              5, false);
  win_map_pub_ =
      nh_.advertise<graph_slam_uk::NDTMapMsg>(win_pub_topic_, 5, false);
  map_pcl_pub_ = nh_.advertise<Pcl>(map_pcl_pub_topic_, 1);
  win_pcl_pub_ = nh_.advertise<Pcl>(win_pcl_pub_topic_, 1);

  const int message_cache = 1;
  if (subscribe_mode_ == "ODOM") {
    laser_sub_.subscribe(nh_, laser_topic_, message_cache);
    odom_sub_.subscribe(nh_, odom_topic_, message_cache);
    // sync messages using approximate algorithm
    odom_sync_.connectInput(odom_sub_, laser_sub_);
    odom_sync_.registerCallback(
        boost::bind(&GraphSlamNode::odom_cb, this, _1, _2));
  } else if (subscribe_mode_ == "POSE") {
    laser_sub_.subscribe(nh_, laser_topic_, message_cache);
    pose_sub_.subscribe(nh_, pose_topic_, message_cache);
    // sync messages using approximate algorithm
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

void GraphSlamNode::start()
{
  ros::spin();
}

void GraphSlamNode::initParameters()
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

  map_pub_topic_ =
      nh_private_.param<std::string>("occupancy_map_topic", "/map");

  map_pcl_pub_topic_ =
      nh_private_.param<std::string>("pcl_map_topic", "/pcl_map");

  win_pub_topic_ =
      nh_private_.param<std::string>("window_map_topic", "/window_map");

  win_pcl_pub_topic_ =
      nh_private_.param<std::string>("window_map_topic_pcl", "/window_pcl");

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

void GraphSlamNode::odom_cb(const nav_msgs::Odometry::ConstPtr &odom,
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

void GraphSlamNode::pose_cb(
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

void GraphSlamNode::laser_cb(const sensor_msgs::LaserScan::ConstPtr &laser)
{
  pose_t new_pose;  // current robot position on the odom_frame
  pcl_ptr_t pcl;    // pcl holding laser scan transformed in base_frame of robot
  ros::Time t = laser->header.stamp;
  if (!prepareAllData(t, laser, new_pose, pcl)) {
    publishTF(t);
    return;
  }
  Eigen::Matrix3d covar;

  covar << 1e-4, 0, 0, 0, 1e-4, 0, 0, 0, 1e-4;
  doAlgorithm(t, new_pose, pcl, covar);
}

bool GraphSlamNode::prepareAllData(const ros::Time &time_stamp,
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
  ROS_DEBUG("GraphSlam: Messages transformed.");
  return true;
}

bool GraphSlamNode::preparePoseData(const ros::Time &time_stamp,
                                    const std::string &laser_frame_id,
                                    pose_t &pose, tf::StampedTransform &tf_base)
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
  // auto odom_trans_mat = eigt::transBtwPoses(last_odom_, pose);
  // if (!movedEnough(odom_trans_mat)) {
  //   ROS_INFO("Graph_slam2d:Robot not moved enough-> droping msgs");
  //   return false;
  // }
  return true;
}

bool GraphSlamNode::prepareLaserData(
    const sensor_msgs::LaserScan::ConstPtr &laser,
    const tf::StampedTransform &tf_base, pcl_ptr_t &pcl)
{
  // project laser message to point cloud class in laser frame_id
  sensor_msgs::PointCloud2 laser_pcl_msg;
  Pcl laser_pcl;
  Pcl laser_pcl_base;
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

void GraphSlamNode::doAlgorithm(const ros::Time &time_stamp, const pose_t &odom,
                                pcl_ptr_t &pcl, const Eigen::Matrix3d &covar)
{
  // prepare initial data in first received msg pair
  if (!is_ready_) {
    is_ready_ = true;
    last_odom_ = odom;
    last_pose_ = pose_t::Zero();
    algorithm_->update(transform_t::Identity(), covar, *pcl, time_stamp);
    publishTF(time_stamp);
    return;
  }
  auto odom_trans = eigt::getTransFromPose(last_odom_).inverse() *
                    eigt::getTransFromPose(odom);
  last_odom_ = odom;

  pose_t current_pose = last_pose_;
  current_pose = algorithm_->update(odom_trans, covar, *pcl, time_stamp);
  last_pose_ = current_pose;
  updateTFTrans(odom, current_pose);

  // PUBLISHING MSGS

  publishTF(time_stamp);
  win_map_pub_.publish(algorithm_->getNDTMap(fixed_frame_));
  occ_map_pub_.publish(algorithm_->getOccupancyGrid(fixed_frame_));
  win_pcl_pub_.publish(algorithm_->getPclMap2(fixed_frame_));
  map_pcl_pub_.publish(algorithm_->getPclMap(fixed_frame_));
  // publish serialized graph
  if (serialize_graph) {
    graph_pub_.publish(algorithm_->getGraphSerialized(fixed_frame_));
    // saveDotGraph();
    ROS_DEBUG("Graph_slam2d:GRAPH Markers published.");
  }
  // prepare data for next iteration
  ++seq_;
}

void GraphSlamNode::publishTF(const ros::Time &time)
{
  // publish map tf transform
  // auto w_trans = eigt::getPoseFromTransform(world_tf_trans_);
  // tf::Transform transform;
  // transform.setOrigin(tf::Vector3(w_trans(0), w_trans(1), 0.0));
  // tf::Quaternion q;
  // q.setRPY(0, 0, w_trans(2));
  // transform.setRotation(q);
  tf_broadcast_.sendTransform(
      tf::StampedTransform(world_tf_trans_, time, fixed_frame_, odom_frame_));
}

void GraphSlamNode::updateTFTrans(const pose_t &odom, const pose_t &slam)
{
  // transform new odometry position with last know SLAM position
  // auto odom_trans = eigt::transformPose(odom, world_tf_trans_);
  // find out how much diverged slam pose and odometric pose in this time
  // update transformation ased on new information
  // std::cout << "!!!!!!!!!!!!!!!!!!!!: "
  //           << eigt::getPoseFromTransform(eigt::transBtwPoses(odom, slam))
  //                  .transpose()
  //           << std::endl;

  // world_tf_trans_ = eigt::transBtwPoses(odom, slam);
  tf::Transform odom_to_base = eigenPoseToTF(odom);
  tf::Transform map_to_base = eigenPoseToTF(slam);
  world_tf_trans_ = tf::Transform(map_to_base * odom_to_base.inverse());
}

Eigen::Matrix3d
GraphSlamNode::arrayToMatrix(const boost::array<double, 36> &array) const
{
  Eigen::Matrix3d covar;
  covar << array[0], array[1], array[5], array[6], array[7], array[11],
      array[30], array[31], array[35];
  return covar;
}

bool GraphSlamNode::movedEnough(const eigt::transform2d_t<double> &trans) const
{
  if (eigt::getAngle(trans) < min_rotation_ &&
      eigt::getDisplacement(trans) < min_displacement_)
    return false;
  return true;
}

tf::Transform GraphSlamNode::eigenPoseToTF(const Eigen::Vector3d &pose) const
{
  geometry_msgs::Pose g_pose;
  g_pose.position.x = pose.x();
  g_pose.position.y = pose.y();

  g_pose.orientation.w = cos(pose.z() * 0.5f);
  g_pose.orientation.z = sin(pose.z() * 0.5f);
  tf::Transform res;
  tf::poseMsgToTF(g_pose, res);
  return res;
}

///////////////// MAIN //////////////////////////////////
int main(int argc, char **argv)
{
  ros::init(argc, argv, "graph_slam");
  ros::NodeHandle n;
  ros::NodeHandle n_private("~");
  NdtSlamAlgortihm algorithm;
  GraphSlamNode slam(n, n_private, algorithm);
  slam.start();

  return 0;
}
