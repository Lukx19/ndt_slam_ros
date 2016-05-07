#ifndef GRAPH_SLAM_UK_GAUSS_NEWTON_OPTIMALIZER2D
#define GRAPH_SLAM_UK_GAUSS_NEWTON_OPTIMALIZER2D

#include <graph_slam_uk/optimizer.h>
#include <graph_slam_uk/graph_slam_interfaces.h>
#include <graph_slam_uk/slam2d_policy.h>
#include <graph_slam_uk/pose_graph.h>
#include <graph_slam_uk/olson_loop_detector.h>
#include <visualization_msgs/MarkerArray.h>
#include <geometry_msgs/Point.h>
#include <std_msgs/ColorRGBA.h>

#include <deque>
#include <dynamic_slam_utils/eigen_tools.h>

#include <iostream>

namespace slamuk
{
template <typename T>
class GaussNewtonOptimalize2d : public IGraphOptimalizer2d<T>
{
  typedef Optimizer<Slam2d_Policy, T> optimizer_t;
  typedef Graph<Slam2d_Policy, T> graph_t;
  typedef OlsonLoopDetector<Slam2d_Policy, T> loop_closurer_t;
  typedef Node<Slam2d_Policy, T> node_t;
  typedef Edge<Slam2d_Policy, T> edge_t;
  typedef Slam2d_Policy Policy;

public:
  GaussNewtonOptimalize2d(IScanmatcher2d &matcher)
    : epsilon_(0.001)
    , iterations_(5)
    , last_node_id_(0)
    , prevlast_node_id_(0)
    , first_node_id_(0)
    , first_node_added_(false)
    , matcher_(&matcher)
    , loop_closure_engine_(loop_closurer_t(graph_, *matcher_))
  {
  }
  // virtual ~GaussNewtonOptimalize2d()
  // {
  // }

  // return true if optimalization changed poses in graph
  virtual bool optimalize();
  virtual bool optimalizeIterationaly();
  virtual double calcTotalGraphError() const;
  virtual size_t addPose(const Eigen::Vector3d &position, T &obj);
  virtual size_t addConstrain(size_t node_id_from, size_t node_id_to,
                              const Eigen::Vector3d &trans,
                              const Eigen::Matrix3d &covar);
  // adds constrain between last two added positions
  virtual size_t addLastConstrain(const Eigen::Vector3d &trans,
                                  const Eigen::Matrix3d &covar);
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
  optimizer_t opt_engine_;
  graph_t graph_;
  IScanmatcher2d *matcher_;
  loop_closurer_t loop_closure_engine_;

  void initializeGrapFromOdom();
  visualization_msgs::MarkerArray createListMarkers(std::string frame_id) const;
  visualization_msgs::MarkerArray createArrowMarkers(std::string frame_id) const;
};

template <typename T>
bool GaussNewtonOptimalize2d<T>::optimalize()
{
  initializeGrapFromOdom();
  return opt_engine_.optimizeGraph(graph_, epsilon_, iterations_);
}

template <typename T>
bool GaussNewtonOptimalize2d<T>::optimalizeIterationaly()
{
  initializeGrapFromOdom();
  return opt_engine_.optimizeGraph(graph_, epsilon_, iterations_);
}

template <typename T>
double GaussNewtonOptimalize2d<T>::calcTotalGraphError() const
{
  return opt_engine_.calcTotalError(graph_);
}

template <typename T>
size_t GaussNewtonOptimalize2d<T>::addPose(const Eigen::Vector3d &position,
                                           T &obj)
{
  size_t id = graph_.addNode(node_t(position, obj));
  if (!first_node_added_) {
    first_node_added_ = true;
    first_node_id_ = id;
    last_node_id_ = id;
  } else {
    prevlast_node_id_ = last_node_id_;
    last_node_id_ = id;
  }
  return id;
}

template <typename T>
size_t GaussNewtonOptimalize2d<T>::addConstrain(size_t node_id_from,
                                                size_t node_id_to,
                                                const Eigen::Vector3d &trans,
                                                const Eigen::Matrix3d &covar)
{
  edge_t e(&graph_.getNode(node_id_from), &graph_.getNode(node_id_to), trans,
           covar);
  size_t id = graph_.addEdge(std::move(e));
  graph_.getEdge(id).setState(edge_t::State::ACTIVE);
  return id;
}

template <typename T>
size_t GaussNewtonOptimalize2d<T>::addLastConstrain(
    const Eigen::Vector3d &trans, const Eigen::Matrix3d &covar)
{
  ROS_INFO_STREAM("Adding constrain betwen nodes:" << prevlast_node_id_ << ":"
                                                   << last_node_id_);
  return addConstrain(prevlast_node_id_, last_node_id_, trans, covar);
}

template <typename T>
bool GaussNewtonOptimalize2d<T>::tryLoopClose()
{
  return loop_closure_engine_.tryLoopClose(last_node_id_);
  // return false;
}

template <typename T>
bool GaussNewtonOptimalize2d<T>::tryLoopClose(size_t node_id)
{
  return loop_closure_engine_.tryLoopClose(node_id);
  // return false;
}

template <typename T>
const Eigen::Vector3d &
GaussNewtonOptimalize2d<T>::getPoseLocation(size_t node_id) const
{
  return graph_.getNode(node_id).getPose();
}

template <typename T>
const T &GaussNewtonOptimalize2d<T>::getPoseData(size_t node_id) const
{
  return graph_.getNode(node_id).getDataObj();
}

template <typename T>
const Eigen::Vector3d &
GaussNewtonOptimalize2d<T>::getConstrainTransform(size_t edge_id) const
{
  return graph_.getEdge(edge_id).getTransform();
}

template <typename T>
const Eigen::Matrix3d &
GaussNewtonOptimalize2d<T>::getConstrainInformMat(size_t edge_id) const
{
  return graph_.getEdge(edge_id).getInformationMatrix();
}

template <typename T>
std::pair<size_t, size_t>
GaussNewtonOptimalize2d<T>::getConstrainPoses(size_t edge_id) const
{
  std::pair<size_t, size_t> out;
  out.first = graph_.getEdge(edge_id).getFrom()->getId();
  out.second = graph_.getEdge(edge_id).getTo()->getId();
  return out;
}

template <typename T>
void GaussNewtonOptimalize2d<T>::setEuclideanMaxError(double epsilon)
{
  epsilon_ = epsilon;
}

template <typename T>
void GaussNewtonOptimalize2d<T>::setMaxIterations(size_t count)
{
  iterations_ = count;
}

template <typename T>
visualization_msgs::MarkerArray
GaussNewtonOptimalize2d<T>::getGraphSerialized(std::string frame_id) const
{
  return createArrowMarkers(frame_id);
}

template <typename T>
void GaussNewtonOptimalize2d<T>::getGraphSerialized(std::ostream &stream) const
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
void GaussNewtonOptimalize2d<T>::initializeGrapFromOdom()
{
  if (!first_node_added_)
    return;

  Policy::Pose pose;
  pose << 0, 0, 0;
  graph_.getNode(first_node_id_).setPose(pose);
  edge_t *next_e = nullptr;
  size_t next_nd;
  for (edge_t *e : graph_.getNode(first_node_id_).getEdgesOut()) {
    if (e->getType() == edge_t::Type::ODOM) {
      next_e = e;
      break;
    }
  }
  while (next_e != nullptr) {
    next_nd = next_e->getTo()->getId();
    graph_.getNode(next_nd)
        .setPose(eigt::transformPose(pose, next_e->getTransMatrix()));
    next_e = nullptr;
    for (edge_t *e : graph_.getNode(next_nd).getEdgesOut()) {
      if (e->getType() == edge_t::Type::ODOM) {
        next_e = e;
        break;
      }
    }
  }
}

template <typename T>
visualization_msgs::MarkerArray
GaussNewtonOptimalize2d<T>::createArrowMarkers(std::string frame_id) const
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
      case edge_t::Type::ODOM:
        colors = red;
        break;
      case edge_t::Type::LOOP:
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
GaussNewtonOptimalize2d<T>::createListMarkers(std::string frame_id) const
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
