#include <nav_msgs/Odometry.h>
#include <ros/ros.h>
#include <tf/tf.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_listener.h>

class OdometryTF
{
public:
  OdometryTF(ros::NodeHandle &nh, ros::NodeHandle &nh_private)
  {
    // find tf_prefix if exists add it to tf_prefix_ variable
    std::string tf_prefix = "";
    std::string tf_prefix_path;
    if (nh_private.searchParam("tf_prefix", tf_prefix_path)) {
      nh_private.getParam(tf_prefix_path, tf_prefix);
    }

    nh_private.param<std::string>("base_frame_id", base_frame_, "base_"
                                                                "footprint");
    nh_private.param<std::string>("fixed_frame_id", fixed_frame_, "odom_comb");
    nh_private.param<std::string>("top_frame_id", top_frame_, "odom");

    if (tf_prefix != "") {
      base_frame_ = tf_prefix + "/" + base_frame_;
      top_frame_ = tf_prefix + "/" + top_frame_;
      fixed_frame_ = tf_prefix + "/" + fixed_frame_;
    }

    odom_sub_ = nh.subscribe<nav_msgs::Odometry>(
        "odom", 1, boost::bind(&OdometryTF::odomCb, this, _1));
  }

private:
  std::string fixed_frame_;
  std::string top_frame_;
  std::string base_frame_;
  ros::Subscriber odom_sub_;
  tf::TransformListener tf_list_;
  tf::TransformBroadcaster tf_broadcast_;

  void odomCb(const nav_msgs::OdometryConstPtr &msg)
  {
    tf::StampedTransform top_base_tf;
    try {
      tf_list_.waitForTransform(top_frame_, base_frame_, msg->header.stamp,
                                ros::Duration(1.0));
      tf_list_.lookupTransform(top_frame_, base_frame_, msg->header.stamp,
                               top_base_tf);
    } catch (...) {
      ROS_ERROR_STREAM("OdometryTFPublisher: Not sucessfull in retrieving tf "
                       "tranform base_frame -> top_frame");
      return;
    }

    tf::Pose robot_pose;
    tf::poseMsgToTF(msg->pose.pose, robot_pose);
    tf::Transform fixed_to_top;
    //   fixed_to_base * odom_to_base.inverse
    fixed_to_top = robot_pose * top_base_tf.inverse();
    tf_broadcast_.sendTransform(tf::StampedTransform(
        fixed_to_top, msg->header.stamp, fixed_frame_, top_frame_));
  }
};
int main(int argc, char **argv)
{
  ros::init(argc, argv, "odom_tf_publisher");
  ros::NodeHandle nh;
  ros::NodeHandle nh_private("~");
  OdometryTF od(nh, nh_private);
  ros::spin();

  return 0;
}
