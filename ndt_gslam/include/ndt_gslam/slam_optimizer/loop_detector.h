#ifndef NDT_GSLAM_LOOP_DETECTOR
#define NDT_GSLAM_LOOP_DETECTOR

#include <ndt_gslam/slam_optimizer/graph_slam_interfaces.h>
#include <ndt_gslam/slam_optimizer/pose_graph.h>
#include <ndt_gslam/utils/eigen_tools.h>
#include <pcl/common/centroid.h>
#include <ros/ros.h>
#include <queue>

namespace slamuk
{
namespace internal
{
struct EdgeCov {
  typedef Eigen::Matrix3d Covar;
  typedef Eigen::Transform<double, 2, Eigen::TransformTraits::Affine> Transform;
  size_t start_id_;
  size_t node_id_;  // id of node to search
  Covar cov_;       // cumulative covariance
  Transform t_;     // cumulative transformation
  double distance_;
  size_t hop_;
  EdgeCov()
    : node_id_(0)
    , cov_(Covar::Zero())
    , t_(Transform::Identity())
    , distance_(0)
    , hop_(0)
    , det_cov_(0)
  {
  }
  explicit EdgeCov(size_t start_id)
    : start_id_(start_id)
    , cov_()
    , t_(Transform::Identity())
    , hop_(0)
    , det_cov_(0)
  {
  }

  bool operator<(const EdgeCov &other) const
  {
    return det_cov_ < other.det_cov_ ? true : false;
  }
  EdgeCov operator+(const EdgeCov &new_e) const
  {
    EdgeCov e;
    e.start_id_ = this->start_id_;
    e.hop_ = this->hop_ + 1;
    e.node_id_ = new_e.node_id_;
    double angle = std::atan2(t_.rotation()(1, 0), t_.rotation()(0, 0));
    double co = std::cos(angle);
    double si = std::sin(angle);
    e.cov_ =
        computeJacc1(co, si, new_e.t_.translation()) * this->cov_ *
            computeJacc1(co, si, new_e.t_.translation()).transpose() +
        computeJacc2(co, si) * new_e.cov_ * computeJacc2(co, si).transpose();
    e.det_cov_ = e.cov_.determinant();
    e.t_ = new_e.t_ * this->t_;
    e.distance_ = this->distance_ + new_e.distance_;
    return e;
  }

private:
  double det_cov_;

  Eigen::Matrix3d computeJacc1(double co, double si,
                               const Eigen::Vector2d &t_j) const
  {
    Eigen::Matrix3d res;
    res << 1, 0, -si * t_j(0) - co * t_j(1), 0, 1, co * t_j(0) - si * t_j(1), 0,
        0, 1;
    return res;
  }
  Eigen::Matrix3d computeJacc2(double co, double si) const
  {
    Eigen::Matrix3d res;
    res << co, -si, 0, si, co, 0, 0, 0, 1;
    return res;
  }
};

struct EdgeCovComparer {
  bool operator()(const EdgeCov &a, const EdgeCov &b)
  {
    return a < b;
  }
};

struct ScanInfo {
  double radius_;
  bool is_ready_;

private:
  char padding[7];

public:
  Eigen::Vector2d centroid_;

  ScanInfo() : radius_(), is_ready_(false), centroid_()
  {
  }

  bool isReady()
  {
    return is_ready_;
  }
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
}  // end of internal namespace

template <typename Policy>
struct LoopClosure {
  std::pair<size_t, size_t> vertices_;
  typename Policy::InformMatrix information_;  // cumulative covariance
  typename Policy::TransformMatrix t_;         // cumulative transformation

  LoopClosure(size_t node_from, size_t node_to,
              typename Policy::InformMatrix inf,
              typename Policy::TransformMatrix trans)
    : vertices_(std::make_pair(node_from, node_to))
    , information_(inf)
    , t_(trans)
  {
  }
};

/////////////////////////////////////////////////LOOP CLOSURE DETECTOR CLASS
/////////////
template <typename P, typename T>
class LoopDetector
{
private:
  typedef typename Graph<P, T>::E Edge;
  typedef typename Graph<P, T>::N Node;
  typedef typename P::Id Id;

public:
  LoopDetector(Graph<P, T> *graph, IScanmatcher2d<T> *matcher)
    : graph_(graph)
    , matcher_(matcher)
    , min_dist_(20)
    , max_dist_(15)
    , max_search_depth_(50)
  {
  }

  /**
   * @brief      generates loop closures
   *
   * @param[in]  node_id  The id of the node for beginign of the search
   *
   * @return     ids of the ew edges
   */
  std::vector<LoopClosure<P>> genLoopClosures(Id node_id);

  /**
   * @brief      Adds loop closures to the pose graph.
   *
   * @param[in]  edges  The edges
   *
   * @return     ids of the ew edges
   */
  std::vector<Id> addToGraph(const std::vector<LoopClosure<P>> &edges);

  /**
   * @brief      Removes a from the pose grap.
   *
   * @param[in]  edges  The edges
   */
  void removeFromGraph(const std::vector<Id> &edges);

  /**
   * @brief      Sets the minimum loop distance.
   *
   * @param[in]  dist  The distance
   */
  void setMinLoopDistance(float dist)
  {
    min_dist_ = dist;
  }
  /**
   * @brief      Sets the maximum loop distance.
   *
   * @param[in]  dist  The distance
   */
  void setMaxLoopDistance(float dist)
  {
    max_dist_ = dist;
  }

protected:
  Graph<P, T> *graph_;
  IScanmatcher2d<T> *matcher_;

  double min_dist_;
  double max_dist_;
  size_t max_search_depth_;
  std::vector<internal::ScanInfo> laser_range_;

  std::vector<internal::EdgeCov>
  genPossibleEdges(Id node_id, const internal::EdgeCov &last_edge);

  const internal::ScanInfo &calcScanParams(Id node_id);
  bool isCloseMatch(const internal::EdgeCov &edge);
};

////////////////////////////////////IMPLEMENTATION
////////////////////////////////////
/** Using Dijkstra projection to distribute covariance throught the graph.
 * Overlapping covariances are considered possivle loop closure match. These
 * hyphotheses are scan matched and resulting transformation is recorded.*/

template <typename P, typename T>
std::vector<LoopClosure<P>> LoopDetector<P, T>::genLoopClosures(Id node_id)
{
  typedef internal::EdgeCov EdgeCov;
  ROS_INFO_STREAM("Finding loop closures for node: " << node_id);
  std::vector<LoopClosure<P>> all_constrains;
  size_t nodes_visited = 0;
  bool far_enough = false;

  std::priority_queue<EdgeCov, std::vector<EdgeCov>, internal::EdgeCovComparer>
      fringe;
  graph_->getNode(node_id).setVisited(true);
  ++nodes_visited;
  // initialize start node edges
  auto first_edges = genPossibleEdges(node_id, EdgeCov(node_id));
  for (EdgeCov &e : first_edges) {
    fringe.push(e);
  }
  size_t depth = 0;
  // repeat while some nodes are unvisited
  while (graph_->nodeCount() > nodes_visited && depth < max_search_depth_) {
    size_t curr_node_id = fringe.top().node_id_;
    std::cout << "node: " << curr_node_id << std::endl;
    Node &target_n = graph_->getNode(node_id);
    Node &source_n = graph_->getNode(curr_node_id);
    // enough overlap
    // calculate l2 distance
    double dist =
        (target_n.getPose().head(2) - source_n.getPose().head(2)).norm();
    std::cout << "distance l2:" << dist << std::endl;
    if (dist < max_dist_ && isCloseMatch(fringe.top())) {
      // increase depth. Prevents too many scanmatches
      ++depth;

      // try to match laser scans
      auto guess_trans = eigt::getTransFromPose(target_n.getPose()).inverse() *
                         eigt::getTransFromPose(source_n.getPose());
      guess_trans(0, 2) = guess_trans(1, 2) = 0;

      MatchResult res = matcher_->match(
          source_n.getDataObj(), target_n.getDataObj(), guess_trans.matrix());
      // test if sucesfull match
      if (res.success_) {
        // create edge
        ROS_INFO_STREAM("loop closure: " << node_id << " -> " << curr_node_id
                                         << "accepted by matching");
        all_constrains.push_back(
            LoopClosure<P>(node_id, curr_node_id, res.inform_, res.transform_));
        std::cout << "source: " << source_n.getPose().transpose() << std::endl;
        std::cout << "target: " << target_n.getPose().transpose() << std::endl;
        std::cout << "calc_trans: "
                  << eigt::getPoseFromTransform(res.transform_).transpose()
                  << std::endl;
        std::cout << "guess: "
                  << eigt::getPoseFromTransform(guess_trans).transpose()
                  << std::endl;

      } else {
        ROS_INFO_STREAM("loop closure: " << node_id << " -> " << curr_node_id
                                         << "rejected by matching");
      }
    }
    graph_->getNode(curr_node_id).setVisited(true);
    ++nodes_visited;
    auto f_top = fringe.top();
    fringe.pop();
    auto next_edges = genPossibleEdges(curr_node_id, f_top);
    std::for_each(next_edges.begin(), next_edges.end(),
                  [&fringe](EdgeCov &e) { fringe.push(e); });
  }

  // put whole graph to base state
  std::for_each(graph_->beginNode(), graph_->endNode(),
                [](Node &nd) { nd.setVisited(false); });
  std::for_each(graph_->beginEdge(), graph_->endEdge(),
                [](Edge &e) { e.setUsed(false); });
  return all_constrains;
}

template <typename P, typename T>
std::vector<typename LoopDetector<P, T>::Id>
LoopDetector<P, T>::addToGraph(const std::vector<LoopClosure<P>> &edges)
{
  std::vector<Id> ids;
  ids.reserve(edges.size());
  for (auto &&e : edges) {
    auto &from_n = graph_->getNode(e.vertices_.first);
    auto &to_n = graph_->getNode(e.vertices_.second);
    Edge new_e(&from_n, &to_n, eigt::getPoseFromTransform(e.t_),
               e.information_);
    Id e_id = graph_->addEdge(std::move(new_e));
    graph_->getEdge(e_id).setType(Edge::Type::LOOP);
    graph_->getEdge(e_id).setState(Edge::State::TEMP);
    ids.push_back(e_id);
  }
  return ids;
}

template <typename P, typename T>
void LoopDetector<P, T>::removeFromGraph(const std::vector<Id> &edges)
{
  for (auto &&id : edges) {
    graph_->removeEdge(id);
  }
}

///////////////////////////////////PROTECTED METHODS ///////////////////////

/** Generates all edges going in graph to the past from current node. Outputs
 * only edges heading to unvisited nodes. Propagates covariance and add it up to
 * previously calculated covariance.
*/
template <typename P, typename T>
std::vector<internal::EdgeCov> LoopDetector<P, T>::genPossibleEdges(
    Id node_id, const internal::EdgeCov &last_edge)
{
  std::vector<internal::EdgeCov> outp;
  // travers graph into the past robot positions
  Node &nd = graph_->getNode(node_id);
  for (Edge *edge_in : nd.getEdgesIn()) {
    if (edge_in->isUsed())
      continue;
    edge_in->setUsed(true);
    if (edge_in->getFrom()->isVisited())
      continue;
    internal::EdgeCov ec;
    ec.node_id_ = edge_in->getFrom()->getId();
    ec.cov_ = edge_in->getInformationMatrix().inverse();
    ec.t_ = edge_in->getTransMatrix().inverse();
    ec.distance_ = ec.t_.translation().norm();
    outp.emplace_back(last_edge + ec);
  }
  return outp;
}

template <typename P, typename T>
const internal::ScanInfo &LoopDetector<P, T>::calcScanParams(Id node_id)
{
  // allocate space if not enough for this node
  if (laser_range_.size() <= node_id)
    laser_range_.insert(laser_range_.end(), node_id - laser_range_.size() + 2,
                        internal::ScanInfo());
  // in case of existance of cached data use them
  if (laser_range_[node_id].isReady()) {
    return laser_range_[node_id];
  }
  internal::ScanInfo res;
  auto data = graph_->getNode(node_id).getDataObj();
  res.radius_ = data.getRadius();
  // compute centroid
  res.centroid_ = data.getCentroid();
  res.is_ready_ = true;
  laser_range_[node_id] = res;

  return laser_range_[node_id];
}

/**  Tests if  start node ----> edgeCov ---> current node id overlapping*/
template <typename P, typename T>
bool LoopDetector<P, T>::isCloseMatch(const internal::EdgeCov &edge)
{
  auto a_params = calcScanParams(edge.start_id_);
  auto b_params = calcScanParams(edge.node_id_);

  Eigen::Vector2d delta_c = edge.t_ * b_params.centroid_ - a_params.centroid_;

  Eigen::Vector2d separation =
      std::max(0.0, delta_c.norm() - b_params.radius_ - a_params.radius_) *
      delta_c.normalized();

  Eigen::Matrix2d icov = edge.cov_.inverse().block(0, 0, 2, 2);
  // double m = (separation.transpose() * icov).dot(separation);
  double dist = delta_c.norm();

  if (edge.distance_ < min_dist_)
    return false;
  else
    return true;
}

}  // end of slamuk namespace
#endif
