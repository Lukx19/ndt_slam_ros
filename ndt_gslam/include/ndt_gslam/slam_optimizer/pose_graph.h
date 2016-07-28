#ifndef NDT_GSLAM_POSE_GRAPH
#define NDT_GSLAM_POSE_GRAPH

#include <Eigen/Eigen>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

namespace slamuk
{
template <typename P, typename T>
class Edge;
template <typename P, typename T>
class Node;
template <typename T>
class ListIterator;
template <typename T>
class ConstListIterator;

template <typename P, typename T>
class Graph
{
public:
  typedef Edge<P, T> E;
  typedef Node<P, T> N;
  typedef size_t Id;
  typedef ListIterator<E> EdgeIterator;
  typedef ConstListIterator<E> ConstEdgeIterator;
  typedef ListIterator<N> NodeIterator;
  typedef ConstListIterator<N> ConstNodeIterator;
  Graph()
  {
  }
  Id addEdge(E&& edge);
  Id addNode(N&& node);
  bool removeEdge(size_t id);
  // bool removeNode(const size_t id);
  E& getEdge(const Id id);
  const E& getEdge(const Id id) const;
  N& getNode(const Id id);
  const N& getNode(const Id id) const;
  NodeIterator beginNode();
  ConstNodeIterator cbeginNode() const;
  NodeIterator endNode();
  ConstNodeIterator cendNode() const;
  EdgeIterator beginEdge();
  ConstEdgeIterator cbeginEdge() const;
  EdgeIterator endEdge();
  ConstEdgeIterator cendEdge() const;
  size_t edgeCount();
  size_t nodeCount();
  void printNodes(std::ostream& out);

private:
  std::vector<std::unique_ptr<E>> edges_;
  std::vector<std::unique_ptr<N>> nodes_;
  std::vector<size_t> id_edge_pool_;
  std::vector<size_t> id_node_pool_;
};

template <typename P, typename T>
typename Graph<P, T>::Id Graph<P, T>::addEdge(E&& edg)
{
  std::unique_ptr<E> e_ptr(new E(std::move(edg)));
  size_t id;
  if (id_edge_pool_.empty()) {
    edges_.push_back(std::move(e_ptr));
    size_t edge_id = edges_.size() - 1;
    id = edge_id;
  } else {
    // previously deleted edges left unallocated possitions. Fill in new edge.
    size_t free_id = id_edge_pool_.back();
    id_edge_pool_.pop_back();
    edges_[free_id] = std::move(e_ptr);
    id = free_id;
  }
  edges_[id]->setId(id);
  nodes_[edges_[id]->getFrom()->getId()]->addEdgeOut(edges_[id].get());
  nodes_[edges_[id]->getTo()->getId()]->addEdgeIn(edges_[id].get());
  return id;
}

template <typename P, typename T>
bool Graph<P, T>::removeEdge(Id edge_id)
{
  if (edge_id > edges_.size())
    return false;
  if (std::find(id_edge_pool_.begin(), id_edge_pool_.end(), edge_id) !=
      id_edge_pool_.end())
    return false;
  id_edge_pool_.push_back(edge_id);
  edges_[edge_id].reset(nullptr);
  return true;
}

template <typename P, typename T>
typename Graph<P, T>::Id Graph<P, T>::addNode(N&& node)
{
  // if(id_node_pool_.empty()){
  std::unique_ptr<N> n_ptr(new N(std::move(node)));
  // std::unique_ptr<N> n_ptr(new N);
  nodes_.push_back(std::move(n_ptr));
  size_t node_id = nodes_.size() - 1;
  nodes_[node_id]->setId(node_id);
  return node_id;
  // }else{
  //   // previously deleted edges left unallocated possitions. Fill in new
  //   edge.
  //   size_t free_id=id_node_pool_.back();
  //   nodes_[free_id]=std::forward<N>(node);
  //   id_node_pool_.pop_back();
  //   nodes_[free_id].setId(free_id);
  //   return free_id;
  // }
}

template <typename P, typename T>
typename Graph<P, T>::E& Graph<P, T>::getEdge(const Id id)
{
  if (id >= edges_.size()) {
    throw std::out_of_range("getEdge id is too big");
  } else if (std::find(id_edge_pool_.begin(), id_edge_pool_.end(), id) !=
             id_edge_pool_.end()) {
    throw std::out_of_range("getEdge edge with 'id' is deleted");
  } else {
    return *edges_[id];
  }
}

template <typename P, typename T>
const typename Graph<P, T>::E& Graph<P, T>::getEdge(const Id id) const
{
  if (id >= edges_.size()) {
    throw std::out_of_range("getEdge id is too big");
  } else if (std::find(id_edge_pool_.begin(), id_edge_pool_.end(), id) !=
             id_edge_pool_.end()) {
    throw std::out_of_range("getEdge edge with 'id' is deleted");
  } else {
    return *edges_[id];
  }
}

template <typename P, typename T>
typename Graph<P, T>::N& Graph<P, T>::getNode(const Id id)
{
  if (id >= nodes_.size()) {
    throw std::out_of_range("getNode id is too big");
  } else {
    return *nodes_[id];
  }
}

template <typename P, typename T>
const typename Graph<P, T>::N& Graph<P, T>::getNode(const Id id) const
{
  if (id >= nodes_.size()) {
    throw std::out_of_range("getNode id is too big");
  } else {
    return *nodes_[id];
  }
}

template <typename P, typename T>
typename Graph<P, T>::NodeIterator Graph<P, T>::beginNode()
{
  return NodeIterator(nodes_.begin(), nodes_.end());
}

template <typename P, typename T>
typename Graph<P, T>::ConstNodeIterator Graph<P, T>::cbeginNode() const
{
  return ConstNodeIterator(nodes_.cbegin(), nodes_.cend());
}
template <typename P, typename T>
typename Graph<P, T>::EdgeIterator Graph<P, T>::beginEdge()
{
  return EdgeIterator(edges_.begin(), edges_.end());
}

template <typename P, typename T>
typename Graph<P, T>::ConstEdgeIterator Graph<P, T>::cbeginEdge() const
{
  return ConstEdgeIterator(edges_.cbegin(), edges_.cend());
}

template <typename P, typename T>
typename Graph<P, T>::NodeIterator Graph<P, T>::endNode()
{
  return NodeIterator(nodes_.end(), nodes_.end());
}

template <typename P, typename T>
typename Graph<P, T>::ConstNodeIterator Graph<P, T>::cendNode() const
{
  return ConstNodeIterator(nodes_.cend(), nodes_.cend());
}

template <typename P, typename T>
typename Graph<P, T>::EdgeIterator Graph<P, T>::endEdge()
{
  return EdgeIterator(edges_.end(), edges_.end());
}

template <typename P, typename T>
typename Graph<P, T>::ConstEdgeIterator Graph<P, T>::cendEdge() const
{
  return ConstEdgeIterator(edges_.cend(), edges_.cend());
}

template <typename P, typename T>
size_t Graph<P, T>::edgeCount()
{
  return edges_.size() - id_edge_pool_.size();
}
template <typename P, typename T>
size_t Graph<P, T>::nodeCount()
{
  return nodes_.size() - id_node_pool_.size();
}

template <typename P, typename T>
void Graph<P, T>::printNodes(std::ostream& out)
{
  for (auto nd : nodes_) {
    nd->printNode(out);
  }
}

//////////////////////////////////Graph iterator ////////////////////
template <typename T>
class ListIterator
{
public:
  typedef ListIterator<T> self_type;
  typedef T value_type;
  typedef T& reference;
  typedef T* pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef int difference_type;
  typedef typename std::vector<std::unique_ptr<T>>::iterator DataIter;
  ListIterator()
  {
  }
  ListIterator(const DataIter& begin, const DataIter& end)
    : current_(begin), end_(end)
  {
  }
  self_type operator++()
  {
    ++current_;
    // in case of empty field in list skip this field and go on
    while (*current_ == nullptr && current_ != end_) {
      ++current_;
    }
    return *this;
  }
  self_type operator++(int)
  {
    ++current_;
    // in case of empty field in list skip this field and go on
    while (current_ != end_ && current_->get() == nullptr) {
      ++current_;
    }
    return *this;
  }
  reference operator*()
  {
    return **current_;
  }
  pointer operator->()
  {
    return current_->get();
  }
  bool operator==(const self_type& rhs) const
  {
    return current_ == rhs.current_;
  }
  bool operator!=(const self_type& rhs)
  {
    return !((*this) == rhs);
  }

private:
  DataIter current_;
  DataIter end_;
};

template <typename T>
class ConstListIterator
{
public:
  typedef ConstListIterator<T> self_type;
  typedef T value_type;
  typedef T& reference;
  typedef T* pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef int difference_type;
  typedef typename std::vector<std::unique_ptr<T>>::const_iterator DataIter;
  ConstListIterator()
  {
  }
  ConstListIterator(const DataIter& begin, const DataIter& end)
    : current_(begin), end_(end)
  {
  }
  self_type operator++()
  {
    ++current_;
    // in case of empty field in list skip this field and go on
    while (*current_ == nullptr && current_ != end_) {
      ++current_;
    }
    return *this;
  }
  self_type operator++(int)
  {
    ++current_;
    // in case of empty field in list skip this field and go on
    while (current_ != end_ && current_->get() == nullptr) {
      ++current_;
    }
    return *this;
  }
  reference operator*()
  {
    return **current_;
  }
  pointer operator->()
  {
    return current_->get();
  }
  bool operator==(const self_type& rhs) const
  {
    return current_ == rhs.current_;
  }
  bool operator!=(const self_type& rhs)
  {
    return !((*this) == rhs);
  }

private:
  DataIter current_;
  DataIter end_;
};

///////////////////////////////////NODE//////////////////////////////////////////////////////////
template <typename P, typename T>
class Node
{
  typedef size_t Id;
  typedef std::unique_ptr<T> obj_ptr_t;
  typedef std::vector<Edge<P, T>*> edge_list_t;

public:
  Node() : id_(0), visited_(false), object_(nullptr)
  {
  }
  Node(const typename P::Pose& pose, const T& obj)
    : id_(0), visited_(false), pose_(pose), object_(obj)
  {
  }
  // Node(const Node<P, T>& other)
  // {
  //   std::cout << "Node copy constructor\n";
  //   clone(other);

  // };
  // Node<P,T> & operator=(const Node<P,T> & other){
  //   std::cout << "Node copy assignemnt\n";
  //   clone(other);
  //   return *this;
  // }
private:
  void clone(const Node<P, T>& other)
  {
    std::cout << "before: ";
    printNode(std::cout);
    if (&other == this)
      return;
    id_ = other.id_;
    visited_ = (other.visited_);
    pose_ = (other.pose_);
    object_ = (other.object_);
    edges_out_ = (other.edges_out_);
    edges_in_ = (other.edges_in_);
    printNode(std::cout);
  }

public:
  const typename P::Pose& getPose() const;
  void setPose(const typename P::Pose& pose);
  void addToPose(typename P::Pose&& pose);
  void addToPose(const typename P::Pose& pose);
  bool addEdgeOut(Edge<P, T>* e);
  bool addEdgeIn(Edge<P, T>* e);
  edge_list_t& getEdgesIn();
  edge_list_t& getEdgesOut();
  //    bool removeEdge(Id id);
  //    bool removeEdge(Edge * e);
  void setId(Id id);
  Id getId() const;
  bool isVisited() const;
  void setVisited(bool visit_state);
  T& getDataObj();
  const T& getDataObj() const;
  void addDataObj(T& obj);
  void printNode(std::ostream& out);

private:
  Id id_;
  bool visited_;
  typename P::Pose pose_;
  T object_;
  edge_list_t edges_out_;
  edge_list_t edges_in_;
  // const float PI = 3.1415927f;

  // float constrainAngle(float x);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <typename P, typename T>
const typename P::Pose& Node<P, T>::getPose() const
{
  return pose_;
}

template <typename P, typename T>
void Node<P, T>::setPose(const typename P::Pose& pose)
{
  pose_ = pose;
}
template <typename P, typename T>
void Node<P, T>::addToPose(typename P::Pose&& pose)
{
  pose_ = P::addPoses(pose_, std::move(pose));
}
template <typename P, typename T>
void Node<P, T>::addToPose(const typename P::Pose& pose)
{
  pose_ = P::addPoses(pose_, pose);
}

template <typename P, typename T>
bool Node<P, T>::addEdgeOut(Edge<P, T>* e)
{
  auto it = std::find(edges_out_.begin(), edges_out_.end(), e);
  if (it == edges_out_.end()) {
    // no existing edge like specified id added .
    edges_out_.push_back(e);
    return true;
  }
  return false;
}

template <typename P, typename T>
bool Node<P, T>::addEdgeIn(Edge<P, T>* e)
{
  auto it = std::find(edges_in_.begin(), edges_in_.end(), e);
  if (it == edges_in_.end()) {
    // no existing edge like specified id added .
    edges_in_.push_back(e);
    return true;
  }
  return false;
}

template <typename P, typename T>
typename Node<P, T>::edge_list_t& Node<P, T>::getEdgesIn()
{
  return edges_in_;
}

template <typename P, typename T>
typename Node<P, T>::edge_list_t& Node<P, T>::getEdgesOut()
{
  return edges_out_;
}

// template<typename T,typename P>
// bool Node<P>::removeEdge(Id id){
//     auto it = std::find(edges_.begin(),edges_.end(),id);
//     if(it != edges_.end()){
//         // edge exist ..can be erased.
//         edges_.erase(it);
//         return true;
//     }
//     return false;
// }

template <typename P, typename T>
void Node<P, T>::setId(Id id)
{
  id_ = id;
}

template <typename P, typename T>
typename Node<P, T>::Id Node<P, T>::getId() const
{
  return id_;
}

template <typename P, typename T>
void Node<P, T>::setVisited(bool state)
{
  visited_ = state;
}

template <typename P, typename T>
bool Node<P, T>::isVisited() const
{
  return visited_;
}

template <typename P, typename T>
T& Node<P, T>::getDataObj()
{
  return object_;
}

template <typename P, typename T>
const T& Node<P, T>::getDataObj() const
{
  return object_;
}

template <typename P, typename T>
void Node<P, T>::addDataObj(T& obj)
{
  object_ = obj_ptr_t(*obj);
}

template <typename P, typename T>
void Node<P, T>::printNode(std::ostream& out)
{
  out << "NODE ID:" << id_ << "POSE:" << pose_.transpose() << std::endl;
}

// template<typename P,typename T>
// float Node<P,T>::constrainAngle(float x){
//     x = fmodf(x + PI,PI*2);
//     if (x < 0)
//         x += 2*PI;
//     return x - PI;
// }

/////////////////////////////////////////////////EDGE///////////////////////////////////////
template <typename P, typename T>
class Edge
{
  typedef size_t Id;

public:
  enum State { REMOVED, TEMP, ACTIVE, CONSTRUCTED };
  enum Type { ODOM, LOOP };

public:
  Edge()
    : id_(0)
    , from_(nullptr)
    , to_(nullptr)
    , state_(State::CONSTRUCTED)
    , type_(Type::ODOM)
    , is_used_(false)
  {
  }
  Edge(Node<P, T>* from, Node<P, T>* to, const typename P::Pose& trans,
       const typename P::InformMatrix& inform)
    : id_(0)
    , from_(from)
    , to_(to)
    , state_(State::CONSTRUCTED)
    , type_(Type::ODOM)
    , is_used_(false)
    , transform_(trans)
    , inform_mat_(inform)
  {
  }
  const Node<P, T>* getFrom() const;
  const Node<P, T>* getTo() const;
  void setId(Id new_id);
  size_t getId() const;
  typename P::ErrorVector getError() const;
  const typename P::InformMatrix& getInformationMatrix() const;
  typename P::JacobianPair getJacobianBlocks() const;
  const typename P::Pose& getTransform() const;
  typename P::TransformMatrix getTransMatrix() const;

  void setState(State new_state);
  State getState() const;
  void setType(Type new_type);
  Type getType() const;
  bool isUsed() const;
  void setUsed(bool state);

private:
  Id id_;
  Node<P, T>* from_;
  Node<P, T>* to_;
  State state_;
  Type type_;
  bool is_used_;
  typename P::Pose transform_;
  typename P::InformMatrix inform_mat_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <typename P, typename T>
const Node<P, T>* Edge<P, T>::getTo() const
{
  return to_;
}

template <typename P, typename T>
const Node<P, T>* Edge<P, T>::getFrom() const
{
  return from_;
}

template <typename P, typename T>
void Edge<P, T>::setId(Id new_id)
{
  id_ = new_id;
}

template <typename P, typename T>
size_t Edge<P, T>::getId() const
{
  return id_;
}

template <typename P, typename T>
typename P::ErrorVector Edge<P, T>::getError() const
{
  return P::calcError(from_->getPose(), to_->getPose(), transform_);
}

template <typename P, typename T>
typename P::JacobianPair Edge<P, T>::getJacobianBlocks() const
{
  typename P::Pose xi = from_->getPose();
  typename P::Pose xj = to_->getPose();
  return P::calcJacobianBlocks(xi, xj, transform_);
}

template <typename P, typename T>
const typename P::InformMatrix& Edge<P, T>::getInformationMatrix() const
{
  return inform_mat_;
}

template <typename P, typename T>
const typename P::Pose& Edge<P, T>::getTransform() const
{
  return transform_;
}

template <typename P, typename T>
typename P::TransformMatrix Edge<P, T>::getTransMatrix() const
{
  return P::vecToTransMat(transform_);
}

template <typename P, typename T>
void Edge<P, T>::setState(State new_state)
{
  state_ = new_state;
}

template <typename P, typename T>
typename Edge<P, T>::State Edge<P, T>::getState() const
{
  return state_;
}

template <typename P, typename T>
void Edge<P, T>::setType(Type new_type)
{
  type_ = new_type;
}

template <typename P, typename T>
typename Edge<P, T>::Type Edge<P, T>::getType() const
{
  return type_;
}

template <typename P, typename T>
void Edge<P, T>::setUsed(bool new_state)
{
  is_used_ = new_state;
}

template <typename P, typename T>
bool Edge<P, T>::isUsed() const
{
  return is_used_;
}

}  // end of namespace slamuk
#endif
