#include <dynamic_slam_utils/eigen_tools.h>

eigt::transform2d_t eigt::transBtwPoses(const pose2d_t &from,
                                           const pose2d_t &to)
{
  transform2d_t t;
  Eigen::Rotation2D<double> rot_from(from(2));
  Eigen::Rotation2D<double> rot_to(to(2));
  t.setIdentity();
  t.matrix().block(0, 0, 2, 2) =
      (rot_from.toRotationMatrix().transpose() * rot_to.toRotationMatrix());
  t.matrix().block(0, 2, 2, 1) =
      rot_from.toRotationMatrix() * (to.head(2) - from.head(2));
  return t;
}


eigt::pose2d_t eigt::transformPose(const pose2d_t &pose,
                                      const transform2d_t &trans)
{
  Eigen::Rotation2D<double> rot(pose(2));
  pose2d_t res_pose;
  pose2d_t inc = getPoseFromTransform(trans);
  res_pose.head(2) = pose.head(2) + rot.toRotationMatrix().transpose() * inc.head(2);
  res_pose(2) = normalizeAngle(pose(2) + inc(2));
  return res_pose;
}


eigt::pose2d_t eigt::getPoseFromTransform(const transform2d_t &trans)
{
  pose2d_t pose;
  Eigen::Vector2d translation = trans.translation();
  pose(0) = translation(0);
  pose(1) = translation(1);
  pose(2) = eigt::getAngle(trans);
  return pose;
}


double eigt::getAngle(const transform2d_t &trans)
{
  return std::atan2(trans.rotation()(1, 0), trans.rotation()(0, 0));
}


double eigt::getDisplacement(const transform2d_t &trans)
{
  return std::sqrt(std::pow(trans.matrix()(0, 2), 2.0) +
                   std::pow(trans.matrix()(1, 2), 2.0));
}


double eigt::getAngleDiffrence(const pose2d_t &from, const pose2d_t &to)
{
  return std::atan2(std::sin(to(2) - from(2)), std::cos(to(2) - from(2)));
}


eigt::transform2d_t
eigt::convertToTransform(const Eigen::Matrix<double, 4, 4> &trans)
{
  transform2d_t new_trans;
  new_trans.setIdentity();
  new_trans.matrix().block(0, 0, 2, 2) = trans.block(0, 0, 2, 2);
  new_trans.matrix().block(0, 2, 2, 1) = trans.block(0, 3, 2, 1);
  return new_trans;
}


Eigen::Matrix<double, 4, 4> eigt::convertFromTransform(const transform2d_t &trans)
{
  Eigen::Matrix<double, 4, 4> new_trans;
  pose2d_t pose = getPoseFromTransform(trans);
  Eigen::AngleAxis<double> init_rotation(pose(2), Eigen::Matrix<double, 3, 1>::UnitZ());
  Eigen::Translation<double, 3> init_translation(pose(0), pose(1), 0);
  return (init_translation * init_rotation).matrix();
}


eigt::transform2d_t eigt::getTransFromPose(const pose2d_t &trans)
{
  transform2d_t trans_mtx;
  Eigen::Rotation2D<double> rot(trans(2));
  Eigen::Translation<double, 2> transl(trans.head(2));
  return transl * rot;
}


double eigt::normalizeAngle(double angle)
{
  return atan2(sin(angle), cos(angle));
  //return angle - 2 * M_PI *std::floor((angle + M_PI) / (2 * M_PI));
}
