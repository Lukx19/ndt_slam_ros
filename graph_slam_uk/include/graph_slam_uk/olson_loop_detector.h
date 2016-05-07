#ifndef GRAPH_SLAM_UK_OLSON_LOOP_DETECTOR
#define GRAPH_SLAM_UK_OLSON_LOOP_DETECTOR

#include <ros/ros.h>
#include <graph_slam_uk/pose_graph.h>
#include <dynamic_slam_utils/eigen_tools.h>
#include <graph_slam_uk/graph_slam_interfaces.h>
#include <queue>
#include <Eigen/Eigenvalues>

namespace slamuk
{
template <typename P, typename T>
class OlsonLoopDetector
{
  struct LoopConstrain
  {
    size_t edge_id_;
    size_t hops_count_;
  };

  struct EdgeGroup
  {
    LoopConstrain first_;
    std::vector<LoopConstrain> edges_;

    bool isValid()
    {
      if (edges_.size > 0)
        return true;
      else
        return false;
    }
  };

  struct EdgeCov
  {
    typedef Eigen::Matrix3d covar_t;
    size_t node_id_;  // id of node to search
    covar_t cov_;
    size_t hop_;
    EdgeCov() : node_id_(0), cov_(covar_t::Zero()), hop_(0)
    {
    }
    EdgeCov(size_t node_id, const covar_t &cov)
      : node_id_(node_id), cov_(cov), hop_(1)
    {
    }

    bool operator<(const EdgeCov &other) const
    {
      return cov_.determinant() < other.cov_.determinant() ? true : false;
    }

    static EdgeCov combineEdges(const EdgeCov &orig_e, const EdgeCov &new_e)
    {
      EdgeCov e;
      e.hop_ = orig_e.hop_ + 1;
      e.cov_ = orig_e.cov_ + new_e.cov_;
      e.node_id_ = new_e.node_id_;
      return e;
    }
  };

  struct EdgeCovComparer
  {
    bool operator()(const EdgeCov &a, const EdgeCov &b)
    {
      return a < b;
    }
  };

  typedef Edge<P, T> edge_t;
  typedef Node<P, T> node_t;

public:
  OlsonLoopDetector(Graph<P, T> &graph, IScanmatcher2d &matcher)
    : graph_(&graph), matcher_(&matcher)
  {
  }

  // return true if there were some changes to ACTIVE edges in graph
  bool tryLoopClose(size_t node_id);

protected:
  double max_laser_range_;
  Graph<P, T> *graph_;
  IScanmatcher2d *matcher_;
  std::deque<EdgeGroup> groups_;
  std::vector<float> laser_range_;

  double MATCH_SCORE = 0.1;
  size_t MAX_LOOP_LENGTH = 5;
  double LOCAL_AMBIGUITY_TRESHOLD = 2;
  size_t MAX_DISTANCE_TO_OVERLAP = 3;

  // return number of created groups
  size_t addToGroups(const std::vector<LoopConstrain> &constrains);
  // removes edges which are not possible to use from graph.
  // return new group with only valid edges
  EdgeGroup optimalizeGroup(const EdgeGroup &group);
  Eigen::VectorXi discretize(const Eigen::VectorXd &vec,
                             const Eigen::MatrixXd &consistency) const;
  Eigen::MatrixXd calcConsistencyMatrix(const EdgeGroup &group) const;
  bool getOdomLinkVals(size_t start_node, size_t end_node,
                       typename P::TransformMatrix *trans,
                       typename P::CovarMatrix *covar) const;
  typename P::CovarMatrix
  combineCovar(const typename P::TransformMatrix &trans_a,
               const typename P::TransformMatrix &trans_b,
               const typename P::CovarMatrix &cov_a,
               const typename P::CovarMatrix &cov_b) const;
  //*********************************
  std::vector<LoopConstrain> generateAllClosures(size_t node_id);
  std::vector<EdgeCov> genPossibleEdges(size_t node_id,
                                        const EdgeCov &last_edge);
  double calcLaserRange(size_t node_id);
  double calcMahalanobisDist(size_t start_node, size_t match_node,
                             const Eigen::Matrix3d &covar);
  bool isGlobalyAmbiguious(size_t start_node, size_t match_node,
                           const EdgeCov &edge) const;
  LoopConstrain addLoopClosureEdge(size_t start_node, size_t match_node,
                                   const MatchResult &res, size_t hops);

  size_t getOdomEdgeCount(size_t start_node, size_t end_node) const;
  std::vector<Edge<P, T> *> getOutOdomEdges(size_t node_id) const;
};

template <typename P, typename T>
bool OlsonLoopDetector<P, T>::tryLoopClose(size_t node_id)
{
  // std::vector<size_t> fully_correct_edges;
  // ambiguity testing
  // size_t group_count = addToGroups(generateAllClosures(node_id));
  // if (group_count > 0) {
  //   // there are some closed groups needing optimalization
  //   for (size_t i = 0; i < group_count; ++i) {
  //     EdgeGroup opt_group = optimalizeGroup(groups_.front());
  //     groups_.pop_front();
  //     std::for_each(opt_group.edges_.begin(), opt_group.edges_.end(),
  //                   [&fully_correct_edges](const LoopConstrain &lc) {
  //       fully_correct_edges.push_back(lc.edge_id_);
  //     });
  //   }
  // }
  // // activate valid loop closure edges for graph backend optimization
  // for (size_t edge_id : fully_correct_edges) {
  //   graph_->getEdge(edge_id).setState(edge_t::State::ACTIVE);
  // }
  // if (fully_correct_edges.size() > 0)
  //   return true;
  // else
  //   return false;
  generateAllClosures(node_id);
  return false;
}

// ****************OPTIMIZATION OF LOOP CLOSURE EDGE GROUPS *************
template <typename P, typename T>
size_t OlsonLoopDetector<P, T>::addToGroups(
    const std::vector<LoopConstrain> &constrains)
{
  size_t created_groups = 0;
  if (constrains.empty())
    return 0;
  if (groups_.empty()) {
    // initialize first group
    EdgeGroup gr1;
    gr1.first_ = constrains[0];
    groups_.push_back(gr1);
  }
  for (const auto &constr : constrains) {
    edge_t &test_edge = graph_->getEdge(constr.edge_id_);
    edge_t &major_edge = graph_->getEdge(groups_.back().first_.edge_id_);
    // find how far is test edge from group major - how many odometry edges away
    size_t mxtx_len = getOdomEdgeCount(major_edge.getFrom()->getId(),
                                       test_edge.getFrom()->getId());
    if (mxtx_len == 0)
      mxtx_len = getOdomEdgeCount(test_edge.getFrom()->getId(),
                                  major_edge.getFrom()->getId());

    size_t myty_len = getOdomEdgeCount(test_edge.getTo()->getId(),
                                       major_edge.getTo()->getId());
    if (myty_len == 0)
      myty_len = getOdomEdgeCount(major_edge.getTo()->getId(),
                                  test_edge.getTo()->getId());
    if (mxtx_len + myty_len > MAX_LOOP_LENGTH) {
      // this node in not suitable for this group
      // make it start node of new group
      EdgeGroup new_gr;
      new_gr.first_ = constr;
      new_gr.edges_.push_back(constr);
      groups_.push_back(new_gr);
      ++created_groups;
    } else {
      // this constrain is in same group as major edge
      groups_.back().edges_.push_back(constr);
    }
  }
  return created_groups;
}

template <typename P, typename T>
typename OlsonLoopDetector<P, T>::EdgeGroup
OlsonLoopDetector<P, T>::optimalizeGroup(const EdgeGroup &group)
{
  EdgeGroup new_group;
  Eigen::MatrixXd consistency_mat = calcConsistencyMatrix(group);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(consistency_mat);
  Eigen::VectorXd eigval = solver.eigenvalues();
  Eigen::VectorXd temp = eigval;
  std::sort(temp.data(), temp.data() + temp.size(),
            [](double a, double b) { return a > b; });

  if (temp(0) / temp(1) < LOCAL_AMBIGUITY_TRESHOLD) {
    // discard hyphotheses from this group.
    // in future possible expansion of set could reuse these edges
    new_group.edges_.clear();

  } else {
    // discretization of edge activating vector
    size_t idx = 0;
    eigval.maxCoeff(&idx);
    Eigen::VectorXi activator =
        discretize(solver.eigenvectors().col(idx), consistency_mat);
    // deleting wrong edges plus creating new group with only good edges
    long i = 0;
    for (auto it = group.edges_.begin(); it != group.edges_.end(); ++it) {
      if (activator(i) == 1) {
        new_group.edges_.push_back(*it);
      } else {
        graph_->removeEdge(it->edge_id_);
      }
      ++i;
    }
    if (new_group.edges_.size() > 0)
      new_group.first_ = new_group.edges_[0];
  }

  return new_group;
}

template <typename P, typename T>
Eigen::VectorXi OlsonLoopDetector<P, T>::discretize(
    const Eigen::VectorXd &vec, const Eigen::MatrixXd &consistency) const
{
  Eigen::VectorXd temp_vec = vec;
  Eigen::VectorXd max_vec = vec;
  double max = 0;
  for (long i = 0; i < vec.size(); ++i) {
    double treshold = vec(i);
    for (long k = 0; k < vec.size(); ++k) {
      if (vec(k) < treshold)
        temp_vec(k) = 0;
      else
        temp_vec(k) = 1;
    }
    double lambda = ((temp_vec.transpose() * consistency).dot(temp_vec)) /
                    temp_vec.transpose().dot(temp_vec);
    if (lambda > max) {
      max = lambda;
      max_vec = temp_vec;
    }
  }
  return max_vec.cast<int>();
}

template <typename P, typename T>
Eigen::MatrixXd
OlsonLoopDetector<P, T>::calcConsistencyMatrix(const EdgeGroup &group) const
{
  size_t edge_count = group.edges_.size();
  Eigen::MatrixXd res_matrix(edge_count, edge_count);
  res_matrix.setZero();
  for (size_t i = 0; i < edge_count; ++i) {
    for (size_t j = i + 1; j < edge_count; ++j) {
      edge_t &hyp_a = graph_->getEdge(group.edges_[i].edge_id_);
      edge_t &hyp_b = graph_->getEdge(group.edges_[j].edge_id_);
      typename P::TransformMatrix T1;
      typename P::CovarMatrix T1_cov;
      if (!getOdomLinkVals(hyp_b.getTo()->getId(), hyp_a.getTo()->getId(), &T1,
                           &T1_cov)) {
        ROS_DEBUG("First odom Link not colculated properly - check missing "
                  "odom edges");
        continue;
      }
      typename P::TransformMatrix T2;
      typename P::CovarMatrix T2_cov;
      if (!getOdomLinkVals(hyp_a.getFrom()->getId(), hyp_b.getFrom()->getId(),
                           &T2, &T2_cov)) {
        ROS_DEBUG("Second odom Link not colculated properly - check missing "
                  "odom edges");
        continue;
      }
      typename P::TransformMatrix T_fin = T1;
      typename P::CovarMatrix T_fin_cov = T1_cov;
      // T1 * hyp A
      T_fin_cov =
          combineCovar(T_fin, hyp_a.getTransMatrix().inverse(), T_fin_cov,
                       hyp_a.getInformationMatrix().inverse());
      T_fin = T1 * hyp_a.getTransMatrix().inverse();
      // (T1 * hyp A) * T2
      T_fin_cov = combineCovar(T_fin, T2, T_fin_cov, T2_cov);
      T_fin = T_fin * T2;
      // ((T1 * hyp A) * T2) * hyp B
      T_fin_cov = combineCovar(T_fin, hyp_b.getTransMatrix(), T_fin_cov,
                               hyp_b.getInformationMatrix().inverse());
      T_fin = T_fin * hyp_b.getTransMatrix();
      typename P::Pose T_fin_short = P::transMatToVec(T_fin);
      res_matrix(i, j) = res_matrix(j, i) =
          std::exp(T_fin_short.dot(T_fin_cov.inverse() * T_fin_short));
    }
  }
  return res_matrix;
}

template <typename P, typename T>
bool
OlsonLoopDetector<P, T>::getOdomLinkVals(size_t start_node, size_t end_node,
                                         typename P::TransformMatrix *trans,
                                         typename P::CovarMatrix *covar) const
{
  // CONSTRAIN: only one odometry edge from node xi to xi+1
  if (start_node == end_node) {
    trans->setIdentity();
    covar->setZero();
    return true;
  }

  edge_t *top_e;
  std::vector<edge_t *> out_edges = getOutOdomEdges(start_node);
  if (out_edges.size() > 0)
    top_e = out_edges[0];
  else
    return false;
  typename P::CovarMatrix T_cov = top_e->getInformationMatrix().inverse();
  typename P::TransformMatrix T_trans = top_e->getTransMatrix();

  while (top_e->getTo()->getId() != end_node) {
    out_edges = getOutOdomEdges(top_e->getTo()->getId());
    if (out_edges.size() > 0)
      top_e = out_edges[0];
    else
      return false;
    T_cov = combineCovar(T_trans, top_e->getTransMatrix(), T_cov,
                         top_e->getInformationMatrix().inverse());
    T_trans = T_trans * top_e->getTransMatrix();
  }
  (*trans) = T_trans;
  (*covar) = T_cov;
  return true;
}

template <typename P, typename T>
typename P::CovarMatrix OlsonLoopDetector<P, T>::combineCovar(
    const typename P::TransformMatrix &trans_a,
    const typename P::TransformMatrix &trans_b,
    const typename P::CovarMatrix &cov_a,
    const typename P::CovarMatrix &cov_b) const
{
  double angle_a = eigt::getAngle(trans_a);
  typename P::Pose delta_b = P::transMatToVec(trans_b);
  double cos_a = std::cos(angle_a);
  double sin_a = std::sin(angle_a);
  double xb = delta_b(0);
  double yb = delta_b(1);
  typename P::JacobianMatrix jaccobian_a, jaccobian_b;
  jaccobian_a << 1, 0, -xb *sin_a - yb *cos_a, 0, 1, -yb *sin_a + xb *cos_a, 0,
      0, 1;
  jaccobian_b << cos_a, -sin_a, 0, sin_a, cos_a, 0, 0, 0, 1;
  return jaccobian_a * cov_a * jaccobian_a.transpose() +
         jaccobian_b * cov_b * jaccobian_b.transpose();
}

//********************************** GENERATING LOOP CLOSURE EDGES********
template <typename P, typename T>
std::vector<typename OlsonLoopDetector<P, T>::LoopConstrain>
OlsonLoopDetector<P, T>::generateAllClosures(size_t node_id)
{
  std::vector<LoopConstrain> all_constrains;
  size_t nodes_visited = 0;
  std::priority_queue<EdgeCov, std::vector<EdgeCov>, EdgeCovComparer> fringe;
  graph_->getNode(node_id).setVisited(true);
  ++nodes_visited;
  // initialize start node edges
  auto first_edges = genPossibleEdges(node_id, EdgeCov());
  for (EdgeCov &e : first_edges) {
    fringe.push(e);
  }
  // repeat while some nodes are unvisited
  while (graph_->nodeCount() != nodes_visited) {
    size_t curr_node_id = fringe.top().node_id_;
    // stop consecutive nodes becoming loop closures
    if (fringe.top().hop_ > 1) {
      // enough overlap
      double dist =
          calcMahalanobisDist(node_id, curr_node_id, fringe.top().cov_);
      if (dist < MAX_DISTANCE_TO_OVERLAP) {
        // calculate transformation with scanmatching
        node_t &source_n = graph_->getNode(curr_node_id);
        node_t &target_n = graph_->getNode(node_id);
        // try to match laser scans
        // MatchResult res =
        //     matcher_->match(source_n.getDataObj(), target_n.getDataObj(),
        //                     eigt::transBtwPoses(target_n.getPose(),
        //                                         source_n.getPose()).matrix());
        // // test if sucesfull match
        // if (res.success_ && res.score_ > MATCH_SCORE) {
        //   // create edge
        //   LoopConstrain lp_constr = addLoopClosureEdge(
        //       node_id, curr_node_id, res, fringe.top().hop_);
        //   all_constrains.push_back(std::move(lp_constr));
        // }
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
                [](Node<P, T> &nd) { nd.setVisited(false); });
  std::for_each(graph_->beginEdge(), graph_->endEdge(),
                [](Edge<P, T> &e) { e.setUsed(false); });
  return all_constrains;
}

template <typename P, typename T>
std::vector<typename OlsonLoopDetector<P, T>::EdgeCov>
OlsonLoopDetector<P, T>::genPossibleEdges(size_t node_id,
                                          const EdgeCov &last_edge)
{
  std::vector<EdgeCov> outp;
  // travers graph into the past robot positionsd
  node_t &nd = graph_->getNode(node_id);
  for (edge_t *edge_in : nd.getEdgesIn()) {
    if (edge_in->isUsed())
      continue;
    edge_in->setUsed(true);
    if (edge_in->getFrom()->isVisited())
      continue;
    EdgeCov ec(edge_in->getFrom()->getId(),
               edge_in->getInformationMatrix().inverse());
    outp.emplace_back(EdgeCov::combineEdges(last_edge, ec));
  }
  // traverse graph to the future if there is any
  // for (edge_t *edge_out : nd->getEdgesOut()) {
  //   if (edge_out->isUsed())
  //     continue;
  //   edge_out->isUsed(true);
  //   if (edge_out->getTo().isVisited())
  //     continue;
  //   EdgeCov ec(edge_out->getTo()->getId(),
  //              edge_out->getInformationMatrix().inverse());
  //   outp.emplace_back(EdgeCov::combineEdges(last_edge, ec));
  // }
  return outp;
}

// specific implementation based on structure of data T
template <typename P, typename T>
double OlsonLoopDetector<P, T>::calcLaserRange(size_t node_id)
{
  // allocate space if not enough for this node
  if (node_id > laser_range_.size())
    laser_range_.insert(laser_range_.end(), node_id - laser_range_.size(),
                        -1.0);
  // in case of existance of cached data use them
  if (laser_range_[node_id] > 0.0)
    return laser_range_[node_id];
  std::vector<float> ranges_horizontal;
  std::vector<float> ranges_vertical;
  auto pcl = graph_->getNode(node_id).getDataObj();
  for (auto &pt : pcl->points) {
    ranges_vertical.push_back(std::abs(pt.data[0] - pcl->sensor_origin_(0)));
    ranges_horizontal.push_back(std::abs(pt.data[1] - pcl->sensor_origin_(1)));
  }
  std::sort(ranges_horizontal.begin(), ranges_horizontal.end());
  std::sort(ranges_vertical.begin(), ranges_vertical.end());

  size_t idx = ranges_vertical.size() * (2 / 3);
  size_t idy = ranges_horizontal.size() * (2 / 3);
  laser_range_[node_id] = std::sqrt(std::pow(ranges_vertical[idx], 2) +
                                    std::pow(ranges_horizontal[idy], 2));
  return laser_range_[node_id];
}

template <typename P, typename T>
double OlsonLoopDetector<P, T>::calcMahalanobisDist(
    size_t start_node, size_t match_node, const Eigen::Matrix3d &covar)
{
  float a_range = calcLaserRange(start_node);
  float b_range = calcLaserRange(match_node);
  Eigen::Vector3d delta_c = (graph_->getNode(match_node).getPose() -
                             graph_->getNode(start_node).getPose());
  Eigen::Vector3d separation =
      std::max(0.0, delta_c.norm() - b_range - a_range) * delta_c.normalized();
  return (separation.transpose() * covar.inverse()).dot(separation);
}

template <typename P, typename T>
bool OlsonLoopDetector<P, T>::isGlobalyAmbiguious(size_t start_node,
                                                  size_t curr_node,
                                                  const EdgeCov &edge) const
{
  return false;
}

template <typename P, typename T>
typename OlsonLoopDetector<P, T>::LoopConstrain
OlsonLoopDetector<P, T>::addLoopClosureEdge(size_t start_node,
                                            size_t match_node,
                                            const MatchResult &res, size_t hops)
{
  node_t &from_n = graph_->getNode(start_node);
  node_t &to_n = graph_->getNode(match_node);
  Edge<P, T> new_e(&from_n, &to_n, eigt::getPoseFromTransform(res.transform_),
                   res.inform_);
  size_t e_id = graph_->addEdge(std::move(new_e));
  graph_->getEdge(e_id).setType(edge_t::Type::LOOP);
  graph_->getEdge(e_id).setState(edge_t::State::TEMP);
  LoopConstrain loop_edge;
  loop_edge.edge_id_ = e_id;
  loop_edge.hops_count_ = hops;

  return loop_edge;
}
template <typename P, typename T>
size_t OlsonLoopDetector<P, T>::getOdomEdgeCount(size_t start_node,
                                                 size_t end_node) const
{
  // CONSTRAIN: only one odometry edge from node xi to xi+1
  if (start_node == end_node)
    return 0;

  edge_t *top_e;
  size_t steps = 1;
  std::vector<edge_t *> out_edges = getOutOdomEdges(start_node);
  if (out_edges.size() > 0)
    top_e = out_edges[0];
  else
    return 0;
  while (top_e->getTo()->getId() != end_node) {
    out_edges = getOutOdomEdges(top_e->getTo()->getId());
    if (out_edges.size() > 0)
      top_e = out_edges[0];
    else
      return 0;
  }
  return steps;
}

template <typename P, typename T>
std::vector<Edge<P, T> *>
OlsonLoopDetector<P, T>::getOutOdomEdges(size_t node_id) const
{
  std::vector<Edge<P, T> *> res;
  auto edges = graph_->getNode(node_id).getEdgesOut();
  for (edge_t *e : edges) {
    if (e->getType() == edge_t::Type::ODOM)
      res.push_back(e);
  }
  return res;
}
}

#endif
