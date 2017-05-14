#include <ndt_gslam/inc_scanmatcher_node.h>

IncScanmatcherNode::IncScanmatcherNode(ros::NodeHandle &n,
                                       ros::NodeHandle &n_private)
  : nh_(n)
  , nh_private_(n_private)
  , last_odom_(Transform::Identity())
  , position_(Transform::Identity())
  , unused_trans_(Transform::Identity())
  , unused_odom_(Transform::Identity())
  , projector_()
  , tf_list_()
  , tf_broadcast_()
  , laser_sub_()
  , odom_pub_()
  , pcl_pub_()
  , occ_pub_()
  , seq_(0)
  , fixed_frame_()
  , robot_base_frame_()
  , odom_frame_()
  , tf_prefix_()
  , laser_topic_()
  , win_radius_(40)
  , cell_size_(0.25)
  , min_rotation_(0.2)
  , min_distance_(0.3)
  , initialized_(false)
  , running_window_()
  , inc_matcher_()
{
  initParameters();
  laser_sub_ = nh_.subscribe<sensor_msgs::LaserScan>(
      laser_topic_, 1, boost::bind(&IncScanmatcherNode::laserCb, this, _1));
  pcl_pub_ = nh_.advertise<Pcl>("win_pcl", 1);
  occ_pub_ = nh_.advertise<nav_msgs::OccupancyGrid>("move_map", 1);
  odom_pub_ = nh_.advertise<nav_msgs::Odometry>("odom", 5);
  inc_matcher_.setCellSize(cell_size_);
  running_window_.reset(new FrameType(cell_size_));
  running_window_->enlarge(-win_radius_, -win_radius_, win_radius_,
                           win_radius_);
}

void IncScanmatcherNode::start()
{
  ros::spin();
}

void IncScanmatcherNode::initParameters()
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

  fixed_frame_ =
      nh_private_.param<std::string>("scanmatcher_farme_id", "odom_scan");

  laser_topic_ = "/scan";

  win_radius_ = static_cast<float>(
      nh_private_.param<double>("scanmatch_window_radius", 40.0));

  cell_size_ = static_cast<float>(nh_private_.param<double>("cell_size", 0.25));

  if (tf_prefix_ != "") {
    robot_base_frame_ = tf_prefix_ + "/" + robot_base_frame_;
    odom_frame_ = tf_prefix_ + "/" + odom_frame_;
    fixed_frame_ = tf_prefix_ + "/" + fixed_frame_;
  }
}

void IncScanmatcherNode::laserCb(const sensor_msgs::LaserScan::ConstPtr &laser)
{
  tf::StampedTransform tf_base;
  Transform odom_position;
  if (!preparePoseData(laser->header.stamp, laser->header.frame_id,
                       odom_position, tf_base)) {
    return;
  }
  Pcl::Ptr pcl = prepareLaserData(laser, tf_base);
  if (pcl->size() < 300) {
    return;
  }

  last_odom_ = odom_position;
  if (!initialized_) {
    last_odom_ = last_odom_.inverse() * odom_position;
  }
  Transform real_movement = calcScanMovement(pcl);
  position_ = position_ * real_movement;

  tf::Transform fixed_to_odom;
  //   fixed_to_base * odom_to_base.inverse
  fixed_to_odom = slamuk::eigenPoseToTF(
      eigt::getPoseFromTransform(position_ * last_odom_.inverse()));

  tf_broadcast_.sendTransform(tf::StampedTransform(
      fixed_to_odom, laser->header.stamp, fixed_frame_, odom_frame_));
  publishOdometry(ros::Time::now());
  publishPcl(ros::Time::now());
  publishOccupancyGrid(ros::Time::now());
  ++seq_;
}

IncScanmatcherNode::Pcl::Ptr IncScanmatcherNode::prepareLaserData(
    const sensor_msgs::LaserScan::ConstPtr &laser,
    const tf::StampedTransform &tf_base)
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
    ROS_ERROR_STREAM("IncScanmatcher: Point cloud transformation not successful"
                     << e.what());
  }
  if (laser_pcl_base.size() < 400) {
    ROS_ERROR_STREAM("IncScanmatcher: Not enough points in laser scan "
                     << laser_pcl_base.size());
  }
  // make copy of pointcloud wrapped in shared ptr
  return laser_pcl_base.makeShared();
}

bool IncScanmatcherNode::preparePoseData(const ros::Time &time_stamp,
                                         const std::string &laser_frame_id,
                                         Transform &odom_position,
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
  odom_position = eigt::convertToTransform(trans_robot.matrix());
  return true;
}

IncScanmatcherNode::Transform
IncScanmatcherNode::calcScanMovement(const Pcl::Ptr &pcl)
{
  if (!initialized_) {
    initialized_ = true;
    // enters first scan as it is to running window
    running_window_->initializeSimple(*pcl);
    ROS_INFO_STREAM("[INC SCANMATCHER NODE]: NDT window initialized!");
    return Transform::Identity();
  }

  Pcl::Ptr pcl_out(new Pcl());
  FrameTypePtr local =
      FrameTypePtr(new FrameType(cell_size_, running_window_->getOrigin()));
  local->initializeSimple(*pcl);
  inc_matcher_.setInputTarget(running_window_);
  inc_matcher_.setInputSource(local->getMeans());
  inc_matcher_.align(*pcl_out,
                     eigt::convertFromTransform(position_cumul_).cast<float>());
  ROS_DEBUG_STREAM("[SLAM ALGORITHM]: incremental scanamtching converged:"
                   << inc_matcher_.hasConverged());
  if (!inc_matcher_.hasConverged()) {
    ROS_INFO_STREAM("[SLAM ALGORITHM]: unsucessful scanmatching-> not moving");
    return Transform::Identity();
  }
  // prepare transformation from successful scan registration
  Transform registration_tf = eigt::convertToTransform<double>(
      inc_matcher_.getFinalTransformation().cast<double>());
  position_cumul_ = registration_tf;
  local->transform(registration_tf);
  // merge in new data to running window
  running_window_->mergeIn(*local, true, true);
  // move running window only horizontaly or verticaly if needed
  position_cumul_ = running_window_->move(registration_tf);

  Transform scan_to_base =
      eigt::getTransFromPose(running_window_->getOrigin()) * position_cumul_;

  // calculate pose increase with respect to last known pose in worls of
  // scanmatching
  Transform increase = position_.inverse() * scan_to_base;
  return increase;
}

void IncScanmatcherNode::publishOdometry(const ros::Time &time)
{
  nav_msgs::Odometry msg;
  msg.header.frame_id = fixed_frame_;
  msg.header.seq = seq_;
  msg.header.stamp = time;
  msg.pose.pose = slamuk::EigenToPoseMsg(
      eigt::getPoseFromTransform(position_).cast<float>());
  msg.child_frame_id = robot_base_frame_;
  odom_pub_.publish(msg);
}

void IncScanmatcherNode::publishPcl(const ros::Time &time)
{
  Pcl::Ptr pcl = running_window_->getMeansTransformed();
  pcl->header.frame_id = fixed_frame_;
  pcl->header.stamp = time.toNSec();
  pcl->header.seq = seq_;
  pcl_pub_.publish(pcl);
}

void IncScanmatcherNode::publishOccupancyGrid(const ros::Time &time) const
{
  auto occ = running_window_->createOccupancyGrid(cell_size_ / 10);
  nav_msgs::OccupancyGrid::Ptr occ_msg(new nav_msgs::OccupancyGrid);
  occ_msg->header.frame_id = fixed_frame_;
  occ_msg->header.seq = seq_;
  occ_msg->header.stamp = time;
  occ_msg->info.height = occ.height_;
  occ_msg->info.width = occ.width_;
  occ_msg->info.resolution = occ.resolution_;
  occ_msg->info.origin = slamuk::EigenToPoseMsg(occ.origin_);
  occ_msg->info.map_load_time = time;
  occ_msg->data = std::move(occ.cells_);
  occ_pub_.publish(occ_msg);
}

bool IncScanmatcherNode::movedEnough(const Transform &trans) const
{
  double rotation = std::abs(eigt::getAngle(trans));
  double translation = std::abs(eigt::getDisplacement(trans));
  if (rotation > min_rotation_ || translation > min_distance_)
    return true;
  return false;
}

///////////////// MAIN //////////////////////////////////
int main(int argc, char **argv)
{
  ros::init(argc, argv, "ndt_inc_odom");
  ros::NodeHandle n;
  ros::NodeHandle n_private("~");
  IncScanmatcherNode inc_odom(n, n_private);
  inc_odom.start();

  return 0;
}
