#ifndef GRAPH_SLAM_UK_PCL_NDT_SCANMATCHING
#define GRAPH_SLAM_UK_PCL_NDT_SCANMATCHING

#include <graph_slam_uk/graph_slam_interfaces.h>
#include <ndt_scanmatching2d/ndt2d.h>
#include <dynamic_slam_utils/eigen_tools.h>

#include <ndt_scanmatching2d/d2d_ndt2d_robust.h>

namespace slamuk
{
class PclNdtScanmatcher : public IScanmatcher2d
{
public:
  PclNdtScanmatcher()
  {
  }
  virtual ~PclNdtScanmatcher()
  {
  }

  virtual MatchResult
  match(const pcl_constptr_t &source, const pcl_constptr_t &target,
        const Eigen::Matrix3d &initial_guess = Eigen::Matrix3d::Identity());
  virtual void setGridStep(double step);
  virtual void setMaxRange(double range);
  virtual void setTransformationEpsilon(double epsilon);

protected:
  pcl::D2DNormalDistributionsTransform2DRobust<pcl::PointXYZ, pcl::PointXYZ>
      pcl_matcher_;
};
}

#endif
