#ifndef GRAPH_SLAM_UK_SLAM2D
#define GRAPH_SLAM_UK_SLAM2D

#include <graph_slam_uk/slam_optimizer/graph_slam_interfaces.h>
#include <graph_slam_uk/slam_optimizer/loop_detector.h>
#include <graph_slam_uk/slam_optimizer/pose_graph.h>
#include <graph_slam_uk/slam_optimizer/rrr_g2o_wrapper.h>
#include <graph_slam_uk/slam_optimizer/rrr_loop_proofer.h>
#include <graph_slam_uk/slam_optimizer/slam2d_policy.h>
#include <graph_slam_uk/utils/eigen_tools.h>

#include <geometry_msgs/Point.h>
#include <std_msgs/ColorRGBA.h>
#include <visualization_msgs/MarkerArray.h>

#include <ros/ros.h>

#include <deque>
#include <map>

#include <iostream>

#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/sparse_optimizer.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/types/slam2d/edge_se2.h>
#include <g2o/types/slam2d/vertex_se2.h>

#include <graph_slam_uk/slam_optimizer/max_mixture/types_g2o_mixture.h>

namespace slamuk
{
template <typename T>
class Slam2D : public IGraphOptimalizer2d<T>
{
  typedef Graph<Slam2d_Policy, T> GraphType;
  typedef LoopDetector<Slam2d_Policy, T> LoopDet;
  typedef Node<Slam2d_Policy, T> NodeType;
  typedef Edge<Slam2d_Policy, T> EdgeType;
  typedef Slam2d_Policy Policy;
  typedef g2o::BlockSolver<g2o::BlockSolverTraits<-1, -1>> SlamBlockSolver;
  typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType>
      SlamLinearSolver;
  typedef g2o::VertexSE2 VertexG2O;
  typedef g2o::SE2 PoseG2O;
  // typedef g2o::EdgeSE2 EdgeG2O;
  // typedef VertexSwitchLinear VertexG2O;
  typedef EdgeSE2Mixture EdgeG2O;
  typedef RRRLoopProofer<RRRG2OWrapper> LoopProofer;

public:
  explicit Slam2D(IScanmatcher2d<T> &matcher)
    : epsilon_(0.001)
    , iterations_(5)
    , last_node_id_(0)
    , prevlast_node_id_(0)
    , first_node_id_(0)
    , first_node_added_(false)
    , graph_()
    , matcher_(&matcher)
    , loop_detector_(&graph_, matcher_)
    , g2o_opt_(new g2o::SparseOptimizer())
    , linear_solver_(new SlamLinearSolver())
    , block_solver_(new SlamBlockSolver(linear_solver_.get()))
    , solver_gauss_(
          new g2o::OptimizationAlgorithmGaussNewton(block_solver_.get()))
    , g2o_rrr_wrapper_(g2o_opt_.get())
    , loop_proofer_(&g2o_rrr_wrapper_)
  {
    linear_solver_->setBlockOrdering(false);
    g2o_opt_->setAlgorithm(solver_gauss_.get());
  }

  // virtual ~Slam2D()
  //  {
  //  }

  // return true if optimalization changed poses in graph
  virtual bool optimalize();
  virtual bool optimalizeIterationaly();
  virtual double calcTotalGraphError() const;
  virtual size_t addPose(const Eigen::Vector3d &position, T &obj);
  virtual size_t addConstrain(size_t node_id_from, size_t node_id_to,
                              const Eigen::Vector3d &trans,
                              const Eigen::Matrix3d &inform_mat);
  // adds constrain between last two added positions
  virtual size_t addLastConstrain(const Eigen::Vector3d &trans,
                                  const Eigen::Matrix3d &inform_mat);
  virtual bool tryLoopClose(size_t node_id);
  // try loop close on the last added pose
  virtual bool tryLoopClose();
  virtual const Eigen::Vector3d &getPoseLocation(size_t node_id) const;
  virtual const T &getPoseData(size_t node_id) const;

  virtual const Eigen::Vector3d &getConstrainTransform(size_t edge_id) const;
  virtual const Eigen::Matrix3d &getConstrainInformMat(size_t edge_id) const;
  virtual std::pair<size_t, size_t> getConstrainPoses(size_t edge_id) const;

  virtual void setEuclideanMaxError(double epsilon);
  virtual void setMaxIterations(size_t count);

  virtual visualization_msgs::MarkerArray
  getGraphSerialized(std::string world_frame_id) const;
  virtual void getGraphSerialized(std::ostream &stream) const;

protected:
  double epsilon_;
  size_t iterations_;
  size_t last_node_id_;
  size_t prevlast_node_id_;
  size_t first_node_id_;
  bool first_node_added_;
  GraphType graph_;
  IScanmatcher2d<T> *matcher_;

  std::unique_ptr<g2o::SparseOptimizer> g2o_opt_;  // needed for solving and
                                                   // graph manipulation
  std::unique_ptr<SlamLinearSolver> linear_solver_;
  std::unique_ptr<SlamBlockSolver> block_solver_;
  std::unique_ptr<g2o::OptimizationAlgorithmGaussNewton> solver_gauss_;

  LoopDet loop_detector_;
  RRRG2OWrapper g2o_rrr_wrapper_;
  LoopProofer loop_proofer_;

  std::map<std::pair<size_t, size_t>, size_t> nodes_to_edge_id_;

  void initializeGrapFromOdom();
  void updatePoseGraph();
  visualization_msgs::MarkerArray createListMarkers(std::string frame_id) const;
  visualization_msgs::MarkerArray createArrowMarkers(std::string frame_id) const;
};

//////////////////////////////IMPLEMENTATION ///////////////////////////////

template <typename T>
bool Slam2D<T>::optimalize()
{
  g2o_opt_->initializeOptimization();
  g2o_opt_->optimize(10);
  updatePoseGraph();
  return true;
}

template <typename T>
bool Slam2D<T>::optimalizeIterationaly()
{
  optimalize();
  return true;
}

template <typename T>
double Slam2D<T>::calcTotalGraphError() const
{
  return 0;
}

template <typename T>
size_t Slam2D<T>::addPose(const Eigen::Vector3d &position, T &obj)
{
  // add vertex to pose_graph
  size_t id = graph_.addNode(NodeType(position, obj));
  // add vertex to g2o
  VertexG2O::EstimateType xytheta(
      PoseG2O(position(0), position(1), position(2)));
  VertexG2O *pose = new VertexG2O();
  pose->setId(id);
  pose->setEstimate(xytheta);
  // initialize ids
  if (!first_node_added_) {
    first_node_added_ = true;
    first_node_id_ = id;
    last_node_id_ = id;
    pose->setFixed(true);
  } else {
    prevlast_node_id_ = last_node_id_;
    last_node_id_ = id;
  }
  g2o_opt_->addVertex(pose);
  return id;
}

template <typename T>
size_t Slam2D<T>::addConstrain(size_t node_id_from, size_t node_id_to,
                               const Eigen::Vector3d &trans,
                               const Eigen::Matrix3d &inform_mat)
{
  ROS_INFO_STREAM("[SLAM2D]: Adding constrain betwen nodes:"
                  << node_id_from << "->" << node_id_to);
  EdgeType e(&graph_.getNode(node_id_from), &graph_.getNode(node_id_to), trans,
             inform_mat);
  size_t id = graph_.addEdge(std::move(e));
  graph_.getEdge(id).setState(EdgeType::State::ACTIVE);
  // adding G2O edgeType
  EdgeG2O *edge_g2o = new EdgeG2O();
  edge_g2o->setMeasurement(PoseG2O(trans(0), trans(1), trans(2)));

  edge_g2o->setInformation(inform_mat);

  edge_g2o->vertices()[0] = g2o_opt_->vertex(node_id_from);
  edge_g2o->vertices()[1] = g2o_opt_->vertex(node_id_to);
  edge_g2o->setId(id);
  g2o_opt_->addEdge(edge_g2o);
  nodes_to_edge_id_[std::make_pair(node_id_from, node_id_to)] = id;
  return id;
}

template <typename T>
size_t Slam2D<T>::addLastConstrain(const Eigen::Vector3d &trans,
                                   const Eigen::Matrix3d &inform_mat)
{
  return addConstrain(prevlast_node_id_, last_node_id_, trans, inform_mat);
}

template <typename T>
bool Slam2D<T>::tryLoopClose()
{
  return tryLoopClose(last_node_id_);
  // return false;
}

template <typename T>
bool Slam2D<T>::tryLoopClose(size_t node_id)
{
  std::cout << "loop closing for nodeType: " << node_id << std::endl;
  std::vector<LoopClosure<Slam2d_Policy>> loops =
      loop_detector_.genLoopClosures(node_id);
  // add to graph
  LoopProofer::VerticesPairSet all_loops;
  LoopProofer::VerticesPairSet good_loops;
  LoopProofer::VerticesPairSet bad_loops;
  for (auto &&constrain : loops) {
    LoopProofer::VerticesPair edgeType;
    edgeType.first = constrain.vertices_.first;
    edgeType.second = constrain.vertices_.second;
    std::cout << "potential loop: " << edgeType.first << "~>" << edgeType.second
              << std::endl;
    all_loops.insert(std::move(edgeType));
    size_t id = addConstrain(edgeType.first, edgeType.second,
                             Policy::transMatToVec(constrain.t_),
                             constrain.information_);
    graph_.getEdge(id).setType(EdgeType::Type::LOOP);
  }
  // proof found edges
  // bool has_added = loop_proofer_.optimizeInc(all_loops,bad_loops,good_loops);
  // if(!has_added)
  //  return false;
  // // remove all bad loop closures
  // for (auto &&id_pair : bad_loops) {
  //   // get id of the removable edgeType
  //   size_t edge_id = nodes_to_edge_id_[id_pair];
  //   // remove from pose graph
  //   graph_.removeEdge(edge_id);
  //   auto all_edges = g2o_opt_->vertex(id_pair.first)->edges();
  //   auto edge_iter = std::find_if(all_edges.begin(), all_edges.end(),
  //                                 [edge_id](g2o::HyperGraph::EdgeType *e) {
  //     if (e->id() == edge_id)
  //       return true;
  //     else
  //       return false;
  //   });
  //   // remove from g2o
  //   if (edge_iter != all_edges.end())
  //     g2o_opt_->removeEdge(*edge_iter);
  //   else
  //     ROS_ERROR_STREAM("Loop closure proofer is rejecting edgeType not used
  //     in "
  //                      "g2o graph");
  // }
  // return true;
  return false;
}

template <typename T>
const Eigen::Vector3d &Slam2D<T>::getPoseLocation(size_t node_id) const
{
  return graph_.getNode(node_id).getPose();
}

template <typename T>
const T &Slam2D<T>::getPoseData(size_t node_id) const
{
  return graph_.getNode(node_id).getDataObj();
}

template <typename T>
const Eigen::Vector3d &Slam2D<T>::getConstrainTransform(size_t edge_id) const
{
  return graph_.getEdge(edge_id).getTransform();
}

template <typename T>
const Eigen::Matrix3d &Slam2D<T>::getConstrainInformMat(size_t edge_id) const
{
  return graph_.getEdge(edge_id).getInformationMatrix();
}

template <typename T>
std::pair<size_t, size_t> Slam2D<T>::getConstrainPoses(size_t edge_id) const
{
  std::pair<size_t, size_t> out;
  out.first = graph_.getEdge(edge_id).getFrom()->getId();
  out.second = graph_.getEdge(edge_id).getTo()->getId();
  return out;
}

template <typename T>
void Slam2D<T>::setEuclideanMaxError(double epsilon)
{
  epsilon_ = epsilon;
}

template <typename T>
void Slam2D<T>::setMaxIterations(size_t count)
{
  iterations_ = count;
}

template <typename T>
visualization_msgs::MarkerArray
Slam2D<T>::getGraphSerialized(std::string frame_id) const
{
  return createArrowMarkers(frame_id);
}

template <typename T>
void Slam2D<T>::getGraphSerialized(std::ostream &stream) const
{
  for (auto it = graph_.cbeginNode(); it != graph_.cendNode(); ++it) {
    double x = it->getPose()(0);
    double y = it->getPose()(1);
    stream << "p" << it->getId() << "[ pose = \"" << x << "," << y << "!\"] \n";
  }

  for (auto it = graph_.cbeginEdge(); it != graph_.cendEdge(); ++it) {
    stream << "p" << it->getFrom()->getId() << "->p" << it->getTo()->getId()
           << std::endl;
  }
}

template <typename T>
void Slam2D<T>::initializeGrapFromOdom()
{
  if (!first_node_added_)
    return;

  Policy::Pose pose;
  pose << 0, 0, 0;
  graph_.getNode(first_node_id_).setPose(pose);
  EdgeType *next_e = nullptr;
  size_t next_nd;
  for (EdgeType *e : graph_.getNode(first_node_id_).getEdgesOut()) {
    if (e->getType() == EdgeType::Type::ODOM) {
      next_e = e;
      break;
    }
  }
  while (next_e != nullptr) {
    next_nd = next_e->getTo()->getId();
    graph_.getNode(next_nd).setPose(
        eigt::transformPose(pose, next_e->getTransMatrix()));
    next_e = nullptr;
    for (EdgeType *e : graph_.getNode(next_nd).getEdgesOut()) {
      if (e->getType() == EdgeType::Type::ODOM) {
        next_e = e;
        break;
      }
    }
  }
}

template <typename T>
void Slam2D<T>::updatePoseGraph()
{
  for (auto it = graph_.beginNode(); it != graph_.endNode(); ++it) {
    Eigen::Vector3d new_pose =
        dynamic_cast<VertexG2O *>(g2o_opt_->vertex(it->getId()))
            ->estimate()
            .toVector();
    it->setPose(new_pose);
    it->getDataObj().updatePosition(new_pose);
  }
}

template <typename T>
visualization_msgs::MarkerArray
Slam2D<T>::createArrowMarkers(std::string frame_id) const
{
  visualization_msgs::MarkerArray markers;
  std::array<float, 3> red{{1, 0, 0}};
  std::array<float, 3> green{{0, 1, 0}};
  // add all edges
  size_t id = 0;
  for (auto it = graph_.cbeginEdge(); it != graph_.cendEdge(); ++it) {
    geometry_msgs::Point start;
    start.x = it->getFrom()->getPose()(0);
    start.y = it->getFrom()->getPose()(1);
    start.z = 0;

    geometry_msgs::Point end;
    end.x = it->getTo()->getPose()(0);
    end.y = it->getTo()->getPose()(1);
    end.z = 0;

    visualization_msgs::Marker marker;
    marker.header.frame_id = frame_id;
    marker.header.stamp = ros::Time();
    marker.ns = "slam_graph";
    marker.id = id;
    marker.type = visualization_msgs::Marker::ARROW;
    marker.action = visualization_msgs::Marker::ADD;
    marker.points.push_back(start);
    marker.points.push_back(end);
    marker.scale.x = 0.1;
    marker.scale.y = 0.3;
    marker.scale.z = 1.0;
    marker.color.a = 1.0;
    std::array<float, 3> colors{{1, 0.8f, 1}};
    switch (it->getType()) {
      case EdgeType::Type::ODOM:
        colors = red;
        break;
      case EdgeType::Type::LOOP:
        colors = green;
        break;
    }
    marker.color.r = colors[0];
    marker.color.g = colors[1];
    marker.color.b = colors[2];
    markers.markers.push_back(std::move(marker));
    ++id;
  }
  return markers;
}

template <typename T>
visualization_msgs::MarkerArray
Slam2D<T>::createListMarkers(std::string frame_id) const
{
  visualization_msgs::MarkerArray markers;
  std_msgs::ColorRGBA red;
  red.r = 1.0f;
  red.g = 0;
  red.b = 0;
  red.a = 1.0f;

  std_msgs::ColorRGBA light_red;
  red.r = 0.4f;
  red.g = 0;
  red.b = 0;
  red.a = 1.0f;

  visualization_msgs::Marker marker;
  marker.header.frame_id = frame_id;
  marker.header.stamp = ros::Time();
  marker.ns = "slam_graph";
  marker.id = 0;
  marker.type = visualization_msgs::Marker::LINE_LIST;
  marker.action = visualization_msgs::Marker::ADD;
  marker.scale.x = 0.3;
  marker.pose.position.x = 0;
  marker.pose.position.y = 0;
  marker.pose.position.z = 0;
  marker.pose.orientation.x = 0;
  marker.pose.orientation.y = 0;
  marker.pose.orientation.z = 0;
  marker.pose.orientation.w = 0;
  marker.color.a = 1.0;
  // add all edges
  for (auto it = graph_.cbeginEdge(); it != graph_.cendEdge(); ++it) {
    geometry_msgs::Point start;
    start.x = it->getFrom()->getPose()(0);
    start.y = it->getFrom()->getPose()(1);
    start.z = 0;

    geometry_msgs::Point end;
    end.x = it->getTo()->getPose()(0);
    end.y = it->getTo()->getPose()(1);
    end.z = 0;
    marker.points.push_back(start);
    marker.points.push_back(end);
    marker.colors.push_back(light_red);
    marker.colors.push_back(red);
  }
  markers.markers.push_back(marker);
  return markers;
}

}  // slamuk namespace

#endif
