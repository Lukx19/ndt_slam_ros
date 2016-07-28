#ifndef NDT_GSLAM_EIGEN_TOOLS
#define NDT_GSLAM_EIGEN_TOOLS

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592
#endif

namespace eigt
{
template <typename T = double>
using transform2d_t = Eigen::Transform<T, 2, Eigen::TransformTraits::Affine>;

template <typename T = double>
using point2d_t = Eigen::Matrix<T, 2, 1>;

template <typename T = double>
using pose2d_t = Eigen::Matrix<T, 3, 1>;
// ************************* DECLARATION *************************************

/// Calculates transformation between inputed poses.
template <typename T = double>
transform2d_t<T> transBtwPoses(const pose2d_t<T> &from, const pose2d_t<T> &to);

// Transforms pose in global frame to its new posistion based on transform
template <typename T = double>
pose2d_t<T> transformPose(const pose2d_t<T> &pose,
                          const transform2d_t<T> &trans);

template <typename T = double>
pose2d_t<T> transformConcat(const pose2d_t<T> &pose,
                            const transform2d_t<T> &trans);

template <typename T = double>
transform2d_t<T> transBtwFrames(const pose2d_t<T> &from, const pose2d_t<T> &to);
// maps from transformation matrix 3x3 to vector
// [delta_x,delta_y, detla_angle]
template <typename T = double>
pose2d_t<T> getPoseFromTransform(const transform2d_t<T> &trans);

// returns angle which is rotational part of transformation
template <typename T = double>
T getAngle(const transform2d_t<T> &trans);

// returns square distance from translation part of transform
template <typename T = double>
T getDisplacement(const transform2d_t<T> &trans);

// returns absolute angular diffrence between poses. Selects shorter angle.
template <typename T = double>
T getAngleDiffrence(const pose2d_t<T> &from, const pose2d_t<T> &to);

// Maps 4x4 transformation matrix to 3x3 transformation matrix
template <typename T = double>
transform2d_t<T> convertToTransform(const Eigen::Matrix<T, 4, 4> &trans);

// MAPS from 3x3 transformation matrix to 4x4 Transformation matrix
template <typename T = double>
Eigen::Matrix<T, 4, 4> convertFromTransform(const transform2d_t<T> &trans);

// Maps transformation encoded in pose vector [delta_x,delta_y,delta_angle] to
// transformation matrix
template <typename T = double>
eigt::transform2d_t<T> getTransFromPose(const pose2d_t<T> &trans);

template <typename T = double>
T normalizeAngle(T angle);

template <typename In, typename Out = In>
Eigen::Matrix<Out, 4, 4> vecToMat3d(const Eigen::Matrix<In, 3, 1> &trans);

template <typename In, typename Out = In>
Eigen::Matrix<Out, 3, 3> vecToMat2d(const Eigen::Matrix<In, 3, 1> &trans);

template <typename In, typename Out = In>
Eigen::Matrix<Out, 3, 1> matToVec2d(const Eigen::Matrix<In, 3, 3> &trans);

template <typename In, typename Out = In>
Eigen::Matrix<Out, 3, 1> matToVec2d(const Eigen::Matrix<In, 4, 4> &trans);
}
// ************************* IMPLEMENTATION****************************
template <typename T>
eigt::transform2d_t<T> eigt::transBtwPoses(const pose2d_t<T> &from,
                                           const pose2d_t<T> &to)
{
  pose2d_t<T> trans_pose;
  trans_pose(0) = to(0) - from(0);
  trans_pose(1) = to(1) - from(1);
  trans_pose(2) = normalizeAngle(to(2) - from(2));
  return eigt::getTransFromPose(trans_pose);
}

template <typename T>
eigt::transform2d_t<T> eigt::transBtwFrames(const pose2d_t<T> &from,
                                            const pose2d_t<T> &to)
{
  // calculate base axes of pose from
  Eigen::Matrix<T, 2, 1> k_base0(1, 0);
  Eigen::Matrix<T, 2, 1> k_base1(0, 1);
  auto trans_from = eigt::getTransFromPose(from);
  Eigen::Matrix<T, 2, 1> from_base0 = trans_from.linear() * k_base0;
  Eigen::Matrix<T, 2, 1> from_base1 = trans_from.linear() * k_base1;

  auto trans_to = eigt::getTransFromPose(to);
  Eigen::Matrix<T, 2, 1> to_base0 = trans_to.linear() * k_base0;
  Eigen::Matrix<T, 2, 1> to_base1 = trans_to.linear() * k_base1;

  transform2d_t<T> T1;
  T1.linear() << from_base0, from_base1;

  transform2d_t<T> T2;
  T2.linear() << to_base0, to_base1;

  transform2d_t<T> T_all;
  T_all.linear() = T1.linear() * T2.linear().inverse();
  T_all.translation() = T2.linear().inverse() * (from.head(2) - to.head(2));
  return T_all;
}

template <typename T>
eigt::pose2d_t<T> eigt::transformPose(const pose2d_t<T> &pose,
                                      const transform2d_t<T> &trans)
{
  pose2d_t<T> res_pose;
  pose2d_t<T> inc = getPoseFromTransform(trans);
  res_pose(0) = pose(0) + inc(0);
  res_pose(1) = pose(1) + inc(1);
  res_pose(2) = normalizeAngle(pose(2) + inc(2));
  return res_pose;
}

template <typename T = double>
eigt::pose2d_t<T> eigt::transformConcat(const pose2d_t<T> &pose,
                                        const transform2d_t<T> &trans)
{
  return eigt::getPoseFromTransform(eigt::getTransFromPose(pose) * trans);
}

template <typename T>
eigt::pose2d_t<T> eigt::getPoseFromTransform(const transform2d_t<T> &trans)
{
  pose2d_t<T> pose;
  auto translation = trans.translation();
  pose(0) = translation(0);
  pose(1) = translation(1);
  pose(2) = eigt::getAngle(trans);
  return pose;
}

template <typename T>
T eigt::getAngle(const transform2d_t<T> &trans)
{
  return std::atan2(trans.rotation()(1, 0), trans.rotation()(0, 0));
}

template <typename T>
T eigt::getDisplacement(const transform2d_t<T> &trans)
{
  return std::sqrt(std::pow(trans.matrix()(0, 2), 2.0) +
                   std::pow(trans.matrix()(1, 2), 2.0));
}

template <typename T>
T eigt::getAngleDiffrence(const pose2d_t<T> &from, const pose2d_t<T> &to)
{
  return std::atan2(std::sin(to(2) - from(2)), std::cos(to(2) - from(2)));
}

template <typename T>
eigt::transform2d_t<T>
eigt::convertToTransform(const Eigen::Matrix<T, 4, 4> &trans)
{
  transform2d_t<T> new_trans;
  new_trans.setIdentity();
  new_trans.matrix().block(0, 0, 2, 2) = trans.block(0, 0, 2, 2);
  new_trans.matrix().block(0, 2, 2, 1) = trans.block(0, 3, 2, 1);
  return new_trans;
}

template <typename T>
Eigen::Matrix<T, 4, 4> eigt::convertFromTransform(const transform2d_t<T> &trans)
{
  Eigen::Matrix<T, 4, 4> new_trans;
  pose2d_t<T> pose = getPoseFromTransform<T>(trans);
  Eigen::AngleAxis<T> init_rotation(pose(2), Eigen::Matrix<T, 3, 1>::UnitZ());
  Eigen::Translation<T, 3> init_translation(pose(0), pose(1), 0);
  return (init_translation * init_rotation).matrix();
}

template <typename T>
eigt::transform2d_t<T> eigt::getTransFromPose(const pose2d_t<T> &trans)
{
  transform2d_t<T> trans_mtx;
  Eigen::Rotation2D<T> rot(trans(2));
  Eigen::Translation<T, 2> transl(trans.head(2));
  return transl * rot;
}

template <typename T>
T eigt::normalizeAngle(T angle)
{
  return atan2(sin(angle), cos(angle));
  // return angle - 2 * M_PI *std::floor((angle + M_PI) / (2 * M_PI));
}

template <typename In, typename Out = In>
Eigen::Matrix<Out, 4, 4> eigt::vecToMat3d(const Eigen::Matrix<In, 3, 1> &trans)
{
  Eigen::Matrix<Out, 4, 4> trans_mat = Eigen::Matrix<Out, 4, 4>::Identity();

  trans_mat.block(0, 0, 3, 3).matrix() =
      Eigen::Matrix<Out, 3, 3>(Eigen::AngleAxis<Out>(
          static_cast<Out>(trans(2)), Eigen::Matrix<Out, 3, 1>::UnitZ()));

  trans_mat.block(0, 3, 3, 1).matrix() = Eigen::Matrix<Out, 3, 1>(
      static_cast<Out>(trans(0)), static_cast<Out>(trans(1)), 0.0);

  return trans_mat;
}

template <typename In, typename Out = In>
Eigen::Matrix<Out, 3, 3> eigt::vecToMat2d(const Eigen::Matrix<In, 3, 1> &trans)
{
  Eigen::Matrix<Out, 3, 3> trans_mat = Eigen::Matrix<Out, 3, 3>::Identity();
  trans_mat.block(0, 0, 2, 2).matrix() =
      Eigen::Rotation2D<Out>(static_cast<Out>(trans(2))).toRotationMatrix();
  trans_mat.block(0, 2, 2, 1).matrix() = Eigen::Matrix<Out, 2, 1>(
      static_cast<Out>(trans(0)), static_cast<Out>(trans(1)));

  return trans_mat;
}

template <typename In, typename Out = In>
Eigen::Matrix<Out, 3, 1> matToVec2d(const Eigen::Matrix<In, 3, 3> &trans)
{
  Eigen::Matrix<Out, 3, 1> vec;
  Eigen::Transform<In, 2, Eigen::Affine, Eigen::ColMajor> trans_mat(trans);
  Eigen::Matrix<In, 2, 1> translation = trans_mat.translation();
  vec << static_cast<Out>(translation(0)), static_cast<Out>(translation(1)),
      static_cast<Out>(
          std::atan2(trans_mat.rotation()(1, 0), trans_mat.rotation()(0, 0)));
  return vec;
}

template <typename In, typename Out = In>
Eigen::Matrix<Out, 3, 1> matToVec2d(const Eigen::Matrix<In, 4, 4> &trans)
{
  Eigen::Matrix<Out, 3, 1> vec;
  Eigen::Transform<In, 3, Eigen::Affine, Eigen::ColMajor> trans_mat(trans);
  Eigen::Matrix<In, 3, 1> translation = trans_mat.translation();
  vec << static_cast<Out>(translation(0)), static_cast<Out>(translation(1)),
      static_cast<Out>(
          std::atan2(trans_mat.rotation()(1, 0), trans_mat.rotation()(0, 0)));
  return vec;
}

#endif
