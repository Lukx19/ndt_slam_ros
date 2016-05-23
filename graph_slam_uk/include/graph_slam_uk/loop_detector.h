#ifndef GRAPH_SLAM_UK_LOOP_DETECTOR
#define GRAPH_SLAM_UK_LOOP_DETECTOR

#include <graph_slam_uk/pose_graph.h>
#include <dynamic_slam_utils/eigen_tools.h>
#include <graph_slam_uk/graph_slam_interfaces.h>
#include <queue>
#include <pcl/common/centroid.h>

namespace slamuk
{

namespace internal
{
    struct EdgeCov
    {
      typedef Eigen::Matrix3d Covar;
      typedef Eigen::Transform<double, 2, Eigen::TransformTraits::Affine> Transform;
      size_t start_id_;
      size_t node_id_;  // id of node to search
      Covar cov_; // cumulative covariance
      Transform t_; // cumulative transformation
      size_t hop_;
      EdgeCov() : node_id_(0), cov_(Covar::Zero()),t_(Transform::Identity()), hop_(0)
      {
      }
      explicit EdgeCov(size_t start_id)
        : start_id_(start_id), cov_(), t_(Transform::Identity()),hop_(0)
      {
      }

      bool operator<(const EdgeCov &other) const
      {
        return cov_.determinant() < other.cov_.determinant() ? true : false;
      }
      EdgeCov operator+(const EdgeCov &new_e) const{
        EdgeCov e;
        e.start_id_ = this->start_id_;
        e.hop_ = this->hop_ + 1;
        e.cov_ = this->cov_ + new_e.cov_;
        e.t_ = new_e.t_ * this->t_;
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

  struct ScanInfo{
    double radius_;
    bool is_ready_;
    Eigen::Vector2d centroid_;

    ScanInfo():radius_(),is_ready_(false),centroid_()
    {}

    bool isReady(){
      return is_ready_;
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };
} // end of internal namespace

template <typename Policy>
struct LoopClosure
{
  typename Policy::Id node_from_;
  typename Policy::Id node_to_;
  typename Policy::InformMatrix information_;  // cumulative covariance
  typename Policy::TransformMatrix t_;         // cumulative transformation

  LoopClosure(typename Policy::Id node_from, typename Policy::Id node_to,
              typename Policy::InformMatrix inf,
              typename Policy::TransformMatrix trans)
    : node_from_(node_from), node_to_(node_to), information_(inf), t_(trans)
  {
  }
};


/////////////////////////////////////////////////LOOP CLOSURE DETECTOR CLASS //////////
template <typename P, typename T>
class LoopDetector
{
private:
    typedef typename Graph<P,T>::E Edge;
    typedef typename Graph<P, T>::N Node;
    typedef typename P::Id Id;
public:
    LoopDetector(Graph<P, T> * graph, IScanmatcher2d * matcher)
     : graph_(graph)
     , matcher_(matcher)
     {}
    std::vector<LoopClosure<P>> genLoopClosures(Id node_id);
    std::vector<Id> addToGraph(const std::vector<LoopClosure<P>> & edges);
    void removeFromGraph(const std::vector<Id> & edges);
protected:
    Graph<P, T> *graph_;
    IScanmatcher2d *matcher_;

    const float MATCH_SCORE  = 0;
    std::vector<internal::ScanInfo> laser_range_;

    std::vector<internal::EdgeCov> genPossibleEdges(
        Id node_id, const internal::EdgeCov &last_edge);

    const internal::ScanInfo & calcScanParams(Id node_id);
    bool isCloseMatch(const internal::EdgeCov &edge);
};

////////////////////////////////////IMPLEMENTATION /////////////////////////////////
/** Using Dijkstra projection to distribute covariance throught the graph.
 * Overlapping covariances are considered possivle loop closure match. These
 * hyphotheses are scan matched and resulting transformation is recorded.*/

template <typename P, typename T>
std::vector<LoopClosure<P>>
LoopDetector<P, T>::genLoopClosures(Id node_id)
{
  typedef internal::EdgeCov EdgeCov;
  std::cout<<"Finding loop closures for node: "<<node_id<<std::endl;
  std::vector<LoopClosure<P>> all_constrains;
  size_t nodes_visited = 0;
  bool far_enough = false;

  std::priority_queue<EdgeCov, std::vector<EdgeCov>, internal::EdgeCovComparer> fringe;
  graph_->getNode(node_id).setVisited(true);
  ++nodes_visited;
  // initialize start node edges
  auto first_edges = genPossibleEdges(node_id, EdgeCov(node_id));
  for (EdgeCov &e : first_edges) {
    fringe.push(e);
  }
  // repeat while some nodes are unvisited
  while (graph_->nodeCount() != nodes_visited) {
    size_t curr_node_id = fringe.top().node_id_;
    std::cout<<"on top of fringe: "<<curr_node_id<<"hops: "<<fringe.top().hop_<<std::endl;
    // enough overlap
    if (isCloseMatch(fringe.top())) {
      if(far_enough){
        // calculate transformation with scanmatching
        Node &source_n = graph_->getNode(curr_node_id);
        Node &target_n = graph_->getNode(node_id);
        // try to match laser scans
        MatchResult res =
             matcher_->match(source_n.getDataObj(), target_n.getDataObj(),
                             eigt::transBtwPoses(target_n.getPose(),
                                                 source_n.getPose()).matrix());
        // test if sucesfull match
        if (res.success_ && res.score_ > MATCH_SCORE) {
          // create edge
          all_constrains.push_back(
              LoopClosure<P>(node_id, curr_node_id, res.inform_, res.transform_));
        }
      }
    }else{
      far_enough = true;
    }
    graph_->getNode(curr_node_id).setVisited(true);
    ++nodes_visited;
    auto f_top = fringe.top();
    fringe.pop();
    auto next_edges = genPossibleEdges(curr_node_id, f_top);
    std::for_each(next_edges.begin(), next_edges.end(),
                  [&fringe](EdgeCov &e)
                  {
                    fringe.push(e);
                  });
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
LoopDetector<P, T>::addToGraph(const std::vector<LoopClosure<P>> & edges){
    std::vector<Id> ids;
    ids.reserve(edges.size());
    for(auto && e:edges){
        auto &from_n = graph_->getNode(e.node_from_);
        auto &to_n = graph_->getNode(e.node_to_);
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
void LoopDetector<P, T>::removeFromGraph(const std::vector<Id> & edges){
    for(auto && id:edges){
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
    outp.emplace_back(last_edge + ec);
  }
  return outp;
}

/** Specific implementation based on structure of data object T. This
 *  implementation expect T to have pcl::PointCloud accesible thru operator*.
*/
template <typename P, typename T>
const internal::ScanInfo & LoopDetector<P, T>::calcScanParams(Id node_id)
{
  // allocate space if not enough for this node
  if (laser_range_.size() <= node_id)
    laser_range_.insert(laser_range_.end(), node_id - laser_range_.size()+2,internal::ScanInfo());
  // in case of existance of cached data use them
  if (laser_range_[node_id].isReady()){
    //std::cout<<"used cached data for node: "<<node_id<<" value"<<laser_range_[node_id].radius_<<std::endl;
    return laser_range_[node_id];
  }
  internal::ScanInfo res;
  std::vector<float> ranges_horizontal;
  std::vector<float> ranges_vertical;
  auto pcl = graph_->getNode(node_id).getDataObj();
  for (auto &pt : pcl->points) {
    ranges_vertical.push_back(std::abs(pt.data[0]));
    ranges_horizontal.push_back(std::abs(pt.data[1]));
  }
  // sort lengths 0,1,....n
  std::sort(ranges_horizontal.begin(), ranges_horizontal.end());
  std::sort(ranges_vertical.begin(), ranges_vertical.end());
  // reject 1/3 of the most distante points
  size_t idy = (ranges_vertical.size() / 3) * 2;
  size_t idx = (ranges_horizontal.size() / 3) * 2;
  res.radius_ = std::max(std::abs(ranges_vertical[idy]),std::abs(ranges_horizontal[idx]));
  // compute centroid
  Eigen::Vector4f centroid;
  pcl::compute3DCentroid(*pcl,centroid);
  res.centroid_ = centroid.head(2).cast<double>();
  res.is_ready_=true;
  laser_range_[node_id] = res;
  //std::cout<<"calc laser range for node:"<<node_id<<" radius:"<<laser_range_[node_id].radius_<<std::endl;

  return laser_range_[node_id];
}

/**  Tests if  start node ----> edgeCov ---> current node id overlapping*/
template <typename P, typename T>
bool LoopDetector<P, T>::isCloseMatch(
    const internal::EdgeCov & edge)
{
   auto a_params = calcScanParams(edge.start_id_);
   auto b_params = calcScanParams(edge.node_id_);
  Eigen::Vector2d delta_c = edge.t_ * b_params.centroid_ - a_params.centroid_;
  // std::cout<<"distance:"<<delta_c.norm()<<std::endl;
  //double delta_angle = std::atan2(edge.t_.rotation()(1, 0), edge.t_.rotation()(0, 0));
  //Eigen::Vector3d delta_c(edge.t_.matrix()(0, 2),edge.t_.matrix()(1, 2),delta_angle);

  Eigen::Vector2d separation =
      std::max(0.0, delta_c.norm() - b_params.radius_ - a_params.radius_) * delta_c.normalized();
  //std::cout<<separation.transpose()<<std::endl;
      Eigen::Matrix2d icov =  edge.cov_.inverse().block(0,0,2,2);
  //double dist = (separation.transpose() * icov).dot(separation);
  double dist = delta_c.norm() - b_params.radius_ - a_params.radius_;
   //std::cout<<"distance_raw:"<<delta_c.norm() - b_params.radius_ - a_params.radius_<<std::endl;
  //std::cout<<"distance:"<<dist<<std::endl;
  if(dist < 0)
    return true;
  else
    return false;
}

} // end of slamuk namespace
#endif
