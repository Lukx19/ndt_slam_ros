#ifndef GRAPH_SLAM_UK_GRAPH_OPTIMALIZER2D_BASE
#define GRAPH_SLAM_UK_GRAPH_OPTIMALIZER2D_BASE
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <boost/shared_ptr.hpp>
#include <visualization_msgs/MarkerArray.h>

namespace slamuk
{
template <typename T>
class IGraphOptimalizer2d
{
  // Implementation must be able to get PointCloud<>::Ptr out of the object T
  // pcl should be accesible by dereferencing operator*
  // pointcloud should be in default coordinate frame of node
  // it is needed for implementation of loop closure with use of IScanmatcher2d
  // interface

  // IDs returned must stay consistant. removing or adding Poses or constrains
  // may not change previously assigned IDs
public:
  virtual ~IGraphOptimalizer2d()
  {
  }
  // return true if optimalization changed poses in graph
  virtual bool optimalize() = 0;
  virtual bool optimalizeIterationaly() = 0;
  virtual double calcTotalGraphError() const = 0;
  virtual size_t addPose(const Eigen::Vector3d &position, T &obj) = 0;
  virtual size_t addConstrain(size_t node_id_from, size_t node_id_to,
                              const Eigen::Vector3d &trans,
                              const Eigen::Matrix3d &covar) = 0;
  // adds constrain between last two added positions
  virtual size_t addLastConstrain(const Eigen::Vector3d &trans,
                                  const Eigen::Matrix3d &covar) = 0;
  // return true if any change to graph were made
  virtual bool tryLoopClose(size_t node_id) = 0;
  // try loop close on last edge added to graph
  virtual bool tryLoopClose() = 0;

  virtual const Eigen::Vector3d &getPoseLocation(size_t node_id) const = 0;
  virtual const T &getPoseData(size_t node_id) const = 0;

  virtual const Eigen::Vector3d &getConstrainTransform(size_t edge_id) const = 0;
  virtual const Eigen::Matrix3d &getConstrainInformMat(size_t edge_id) const = 0;
  virtual std::pair<size_t, size_t> getConstrainPoses(size_t edge_id) const = 0;

  virtual void setEuclideanMaxError(double epsilon) = 0;
  virtual void setMaxIterations(size_t count) = 0;

  virtual visualization_msgs::MarkerArray
  getGraphSerialized(std::string world_frame_id) const = 0;
  virtual void getGraphSerialized(std::ostream &stream) const = 0;
};

struct MatchResult
{
public:
  typedef Eigen::Transform<double, 2, Eigen::TransformTraits::Affine>
      transform_t;
  bool success_;
  double score_;
  Eigen::Matrix3d inform_;
  transform_t transform_;

  MatchResult() : success_(false), score_(0.0)
  {
    transform_.setIdentity();
  }

  MatchResult(bool success, double score, const transform_t &trans)
    : success_(success), score_(score), transform_(trans)
  {
  }
};

class IScanmatcher2d
{
public:
  typedef pcl::PointCloud<pcl::PointXYZ> pcl_t;
  typedef pcl_t::Ptr pcl_ptr_t;
  typedef pcl_t::ConstPtr pcl_constptr_t;
  virtual ~IScanmatcher2d()
  {
  }

  virtual MatchResult
  match(const pcl_constptr_t &source, const pcl_constptr_t &target,
        const Eigen::Matrix3d &initial_guess = Eigen::Matrix3d::Identity()) = 0;
  virtual void setGridStep(double step) = 0;
  virtual void setMaxRange(double range) = 0;
  virtual void setTransformationEpsilon(double epsilon) = 0;
};

}  // slamuk namespace
#endif
