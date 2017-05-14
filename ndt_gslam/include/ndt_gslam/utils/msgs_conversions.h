#ifndef NDT_GSLAM_MSGS_CONVERSIONS
#define NDT_GSLAM_MSGS_CONVERSIONS

#include <geometry_msgs/Pose.h>
#include <Eigen/Dense>
namespace slamuk
{
geometry_msgs::Pose EigenToPoseMsg(const Eigen::Vector3f &pose);

geometry_msgs::Pose EigenToPoseMsg(const Eigen::Vector3f &pose)
{
  geometry_msgs::Pose msg;
  msg.position.x = pose(0);
  msg.position.y = pose(1);
  msg.position.z = 0;
  Eigen::Quaternionf quat(Eigen::AngleAxisf(pose(2), Eigen::Vector3f::UnitZ()));
  msg.orientation.x = quat.x();
  msg.orientation.y = quat.y();
  msg.orientation.z = quat.z();
  msg.orientation.w = quat.w();
  return msg;
}
}

#endif
