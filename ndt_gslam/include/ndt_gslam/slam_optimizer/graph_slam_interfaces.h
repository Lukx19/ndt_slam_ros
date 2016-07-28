#ifndef NDT_GSLAM_GRAPH_OPTIMALIZER2D_BASE
#define NDT_GSLAM_GRAPH_OPTIMALIZER2D_BASE
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <visualization_msgs/MarkerArray.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <boost/shared_ptr.hpp>
#include <memory>

namespace slamuk
{
/**
 * @brief      Base interface for graph slam optimization
 *
 * @tparam     T     Class representing measurements data. Must provide methods
 *                   Eigen::Vector2d getCentroid(); const double getRadius();
 *                   const const Data &getData(); const void
 * updatePosition(const
 *                   Eigen::Vector3d &new_pose);
 */
template <class T>
class IGraphOptimalizer2d
{
public:
  virtual ~IGraphOptimalizer2d()
  {
    std::cout << "destructing IGraphOpt" << std::endl;
  }

  /**
   * @brief      Bulk optimization of the whole graph
   *
   * @return     true if the graph has changed
   */
  virtual bool optimalize() = 0;

  /**
   * @brief      Incremental optimization
   *
   * @return     true if the graph has changed
   */
  virtual bool optimalizeIterationaly() = 0;

  /**
   * @brief      Calculates the total graph error.
   *
   * @return     The total graph error.
   */
  virtual double calcTotalGraphError() const = 0;

  /**
   * @brief      Adds a pose to the graph
   *
   * @param[in]  position  The position
   * @param      obj       The measurement object
   *
   * @return     unique id of the node
   */
  virtual size_t addPose(const Eigen::Vector3d &position, T &obj) = 0;

  /**
   * @brief      Adds a constrain.
   *
   * @param[in]  node_id_from  The node identifier from
   * @param[in]  node_id_to    The node identifier to
   * @param[in]  trans         The transformation
   * @param[in]  inform_mat    The information matrix
   *
   * @return     unique id of the constrain
   */
  virtual size_t addConstrain(size_t node_id_from, size_t node_id_to,
                              const Eigen::Vector3d &trans,
                              const Eigen::Matrix3d &inform_mat) = 0;

  /**
   * @brief      Adds a constrain between last two nodes.
   *
   * @param[in]  trans       The transformation
   * @param[in]  inform_mat  The information matrix
   *
   * @return     unique id of the constrain
   */
  virtual size_t addLastConstrain(const Eigen::Vector3d &trans,
                                  const Eigen::Matrix3d &inform_mat) = 0;
  // return true if any change to graph were made

  /**
   * @brief      Loop closure detection and validation
   *
   * @param[in]  node_id  The node identifier where to start
   *
   * @return     true if loop closure added
   */
  virtual bool tryLoopClose(size_t node_id) = 0;

  /**
   * @brief      Loop closure detection and validation starting at last pose in
   * the graph
   *
   * @return     true if loop closure added
   */
  virtual bool tryLoopClose() = 0;

  /**
   * @brief      Gets the pose location.
   *
   * @param[in]  node_id  The node identifier
   *
   * @return     The pose
   */
  virtual const Eigen::Vector3d &getPoseLocation(size_t node_id) const = 0;

  /**
   * @brief      Gets the pose data.
   *
   * @param[in]  node_id  The node identifier
   *
   * @return     The pose data.
   */
  virtual const T &getPoseData(size_t node_id) const = 0;

  /**
   * @brief      Gets the constrain transform.
   *
   * @param[in]  edge_id  The edge identifier
   *
   * @return     The constrain transform.
   */
  virtual const Eigen::Vector3d &getConstrainTransform(size_t edge_id) const = 0;

  /**
   * @brief      Gets the constrain information matrix.
   *
   * @param[in]  edge_id  The edge identifier
   *
   * @return     The constrain information matrix.
   */
  virtual const Eigen::Matrix3d &getConstrainInformMat(size_t edge_id) const = 0;

  /**
   * @brief      Gets the constrain nodes
   *
   * @param[in]  edge_id  The edge identifier
   *
   * @return     <The node id from , node id to>
   */
  virtual std::pair<size_t, size_t> getConstrainPoses(size_t edge_id) const = 0;

  /**
   * @brief      Sets the euclidean maximum error.
   *
   * @param[in]  epsilon  The epsilon
   */
  virtual void setEuclideanMaxError(double epsilon) = 0;

  /**
   * @brief      Sets the maximum iterations.
   *
   * @param[in]  count  The count
   */
  virtual void setMaxIterations(size_t count) = 0;

  /**
   * @brief      Sets the loop generation minimum distance.
   *
   * @param[in]  dist  The distance
   */
  virtual void setLoopGenerationMinDist(float dist) = 0;

  /**
   * @brief      Sets the loop generation maximum distance.
   *
   * @param[in]  dist  The distance
   */
  virtual void setLoopGenerationMaxDist(float dist) = 0;

  /**
   * @brief      Sets the loop registration score.
   *
   * @param[in]  score  The score
   */
  virtual void setLoopRegistrationScore(float score) = 0;

  /**
   * @brief      Gets the graph serialized.
   *
   * @param[in]  world_frame_id  The world frame identifier
   *
   * @return     The graph msg.
   */
  virtual visualization_msgs::MarkerArray
  getGraphSerialized(std::string world_frame_id) const = 0;
  virtual void getGraphSerialized(std::ostream &stream) const = 0;
};

struct MatchResult {
public:
  typedef Eigen::Transform<double, 2, Eigen::TransformTraits::Affine>
      transform_t;
  Eigen::Matrix3d inform_;
  transform_t transform_;
  bool success_;
  double score_;

  MatchResult() : success_(false), score_(0.0)
  {
    transform_.setIdentity();
  }

  MatchResult(bool success, double score, const transform_t &trans)
    : transform_(trans), success_(success), score_(score)
  {
  }
};

template <typename FrameType>
class IScanmatcher2d
{
public:
  typedef FrameType CloudFrame;
  virtual ~IScanmatcher2d()
  {
  }

  /**
   * @brief      Calculate registration
   *
   * @param[in]  source         The source
   * @param[in]  target         The target
   * @param[in]  initial_guess  The initial guess
   *
   * @return     { description_of_the_return_value }
   */
  virtual MatchResult
  match(const CloudFrame &source, const CloudFrame &target,
        const Eigen::Matrix3d &initial_guess = Eigen::Matrix3d::Identity()) = 0;

  /**
   * @brief      Sets the score threshold.
   *
   * @param[in]  score  The score
   */
  virtual void setScoreThreshold(float score);
};

}  // slamuk namespace
#endif
