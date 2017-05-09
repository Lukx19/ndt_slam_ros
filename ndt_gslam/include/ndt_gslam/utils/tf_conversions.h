#ifndef NDT_GSLAM_UTILS_TF_CONVERSIONS_H
#define NDT_GSLAM_UTILS_TF_CONVERSIONS_H

#include <tf/tf.h>

namespace slamuk
{
tf::Transform eigenPoseToTF(const Eigen::Vector3d &pose);

tf::Transform eigenPoseToTF(const Eigen::Vector3d &pose)
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
}
#endif  // NDT_GSLAM_UTILS_TF_CONVERSIONS_H
