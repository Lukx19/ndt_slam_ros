#ifndef GRAPH_SLAM_UK_NDT_SCANMATCHING
#define GRAPH_SLAM_UK_NDT_SCANMATCHING

#include <graph_slam_uk/slam_optimizer/graph_slam_interfaces.h>
#include <graph_slam_uk/utils/eigen_tools.h>
#include <graph_slam_uk/registration/d2d_ndt2d_robust.h>

namespace slamuk
{
template <typename FrameType>
class NdtScanmatcher : public IScanmatcher2d<FrameType>
{
public:
  NdtScanmatcher()
  {
  }
  virtual ~NdtScanmatcher()
  {
  }

  virtual MatchResult
  match(const FrameType &source, const FrameType &target,
        const Eigen::Matrix3d &initial_guess = Eigen::Matrix3d::Identity());

protected:
  pcl::D2DNormalDistributionsTransform2DRobust<pcl::PointXYZ, pcl::PointXYZ>
      matcher;
};

template <typename FrameType>
MatchResult
NdtScanmatcher<FrameType>::match(const FrameType &source,
                                 const FrameType &target,
                                 const Eigen::Matrix3d &initial_guess)
{
  matcher.setInputSource(source);
  matcher.setInputTarget(target);
  // Set initial alignment estimate found using robot odometry.
  Eigen::Matrix<double, 4, 4> init_guess =
      eigt::convertFromTransform(eigt::transform2d_t<double>(initial_guess));
  // Calculating required rigid transform to align the input cloud to the target
  // cloud.
  pcl::PointCloud<pcl::PointXYZ> pcl_out;
  matcher.align(pcl_out, init_guess.cast<float>());

  ROS_INFO_STREAM("PCL_NDT2D:Normal Distributions Transform has converged:"
                  << matcher.hasConverged());
  MatchResult res;
  res.success_ = matcher.hasConverged();
  res.score_ = 1;
  res.inform_ = matcher.getInformMatrix();
  Eigen::Matrix4f trans = matcher.getFinalTransformation();
  res.transform_ = eigt::convertToTransform<double>(trans.cast<double>());
  return res;
}
}

#endif
