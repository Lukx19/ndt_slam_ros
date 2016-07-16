#ifndef GRAPH_SLAM_UK_NDT_SCANMATCHING
#define GRAPH_SLAM_UK_NDT_SCANMATCHING

#include <graph_slam_uk/registration/d2d_ndt2d_robust.h>
#include <graph_slam_uk/slam_optimizer/graph_slam_interfaces.h>
#include <graph_slam_uk/utils/eigen_tools.h>

#include <graph_slam_uk/utils/point_cloud_tools.h>

namespace slamuk
{
template <typename FrameType>
class NdtScanmatcher : public IScanmatcher2d<FrameType>
{
public:
  NdtScanmatcher()
  {
    matcher.setStepSize(0.1);
    matcher.setOulierRatio(0.55);
  }
  virtual ~NdtScanmatcher()
  {
  }

  virtual MatchResult
  match(const FrameType &source, const FrameType &target,
        const Eigen::Matrix3d &initial_guess = Eigen::Matrix3d::Identity());

protected:
  pcl::D2DNormalDistributionsTransform2DRobust<typename FrameType::Point,
                                               typename FrameType::Point>
      matcher;
};

template <typename FrameType>
MatchResult
NdtScanmatcher<FrameType>::match(const FrameType &source,
                                 const FrameType &target,
                                 const Eigen::Matrix3d &initial_guess)
{
  matcher.setInputSource(source.getData());
  matcher.setInputTarget(target.getData());
  // Set initial alignment estimate found using robot odometry.
  Eigen::Matrix<float, 4, 4> init_guess =
      eigt::convertFromTransform(eigt::transform2d_t<double>(initial_guess))
          .cast<float>();
  // Calculating required rigid transform to align the input cloud to the target
  // cloud.
  typename pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_out(
      new pcl::PointCloud<pcl::PointXYZ>());

  matcher.align(*pcl_out, init_guess);
  // pcl::visualizePcl<typename FrameType::Point>(target.getData()->getMeans(),
  //                                              pcl_out);

  ROS_INFO_STREAM(
      "[NDT_SCANMATCHER]:Normal Distributions Transform has converged:"
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
