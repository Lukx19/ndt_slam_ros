#ifndef NDT_GSLAM_NDT_SCANMATCHING
#define NDT_GSLAM_NDT_SCANMATCHING

#include <ndt_gslam/registration/d2d_ndt2d_robust.h>
#include <ndt_gslam/slam_optimizer/graph_slam_interfaces.h>
#include <ndt_gslam/utils/eigen_tools.h>

#include <ndt_gslam/utils/point_cloud_tools.h>
#include <pcl/io/pcd_io.h>
#include <ctime>
#include <sstream>

namespace slamuk
{
template <typename FrameType>
class NdtScanmatcher : public IScanmatcher2d<FrameType>
{
public:
  NdtScanmatcher() : seq_(0), time_((unsigned)time(0))
  {
    matcher.setStepSize(0.1);
    matcher.setOulierRatio(0.55);
    matcher.setRejectionLimit(0.6);
    matcher.setRotationRange(3);
  }
  virtual ~NdtScanmatcher()
  {
  }

  virtual MatchResult
  match(const FrameType &source, const FrameType &target,
        const Eigen::Matrix3d &initial_guess = Eigen::Matrix3d::Identity());
  virtual void setScoreThreshold(float score) override
  {
    matcher.setRejectionLimit(score);
  }

private:
  size_t seq_;
  size_t time_;
  float score_threshold_;
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
  // Eigen::Matrix<float, 4, 4> init_guess = Eigen::Matrix4f::Identity();
  Eigen::Matrix<float, 4, 4> init_guess =
      eigt::convertFromTransform(eigt::transform2d_t<double>(initial_guess))
          .cast<float>();
  // Calculating required rigid transform to align the input cloud to the target
  // cloud.
  typename pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_out(
      new pcl::PointCloud<pcl::PointXYZ>());

  matcher.align(*pcl_out, init_guess);
  std::string filename;
  std::stringstream ss;
  if (matcher.hasConverged()) {
    ss << "valid_";
  } else {
    ss << "invalid_";
  }
  // save measured point clouds and solutions to the disk for testing
  ss << time_ << "_" << seq_ << "score:" << matcher.getAlignmentQuality();
  ++seq_;
  pcl::savePcl<typename FrameType::Point>(target.getData()->getMeans(), pcl_out,
                                          ss.str());
  std::cout << ss.str() << std::endl;
  std::stringstream ss_source;
  std::stringstream ss_target;

  ss_target << "target_" << time_ << seq_ << ".pcd";
  ss_source << "source_" << time_ << seq_ << ".pcd";
  pcl::io::savePCDFile(ss_target.str(), *(target.getData()->getMeans()));
  pcl::io::savePCDFile(ss_source.str(), *(source.getData()->getMeans()));

  ROS_INFO_STREAM(
      "[NDT_SCANMATCHER]:Normal Distributions Transform has converged:"
      << matcher.hasConverged());
  MatchResult res;
  res.success_ = matcher.hasConverged();
  res.score_ = matcher.getAlignmentQuality();
  res.inform_ =
      (Eigen::Matrix3d::Identity() * 5000 /* matcher.getAlignmentQuality()*/);
  res.inform_(2, 2) = 10000; /* matcher.getAlignmentQuality();*/
  // matcher.getInformMatrix();
  Eigen::Matrix4f trans = matcher.getFinalTransformation();
  res.transform_ = eigt::convertToTransform<double>(trans.cast<double>());
  return res;
}
}

#endif
