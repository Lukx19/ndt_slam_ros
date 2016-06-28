#ifndef GRAPH_SLAM_UK_RRR_G2O_WRAPPER
#define GRAPH_SLAM_UK_RRR_G2O_WRAPPER

#include <g2o/core/sparse_optimizer.h>

namespace slamuk{

//template<typename VertexType,typename EdgeType>
class RRRG2OWrapper{
public:
  typedef std::pair<int, int> EdgeId;
  typedef std::set<EdgeId> EdgeIdSet;
  typedef std::map<std::pair<int, int>, double> LoopErrorMap;

private:
  typedef g2o::SparseOptimizer::Edge* EdgePtr;
  typedef g2o::SparseOptimizer::EdgeSet EdgePtrSet;
  typedef std::map<EdgeId, EdgePtr> EdgeIdToEdgePtrMap;

public:
  struct Result
  {
    Result() : edge_count_(0), graph_error_(-1)
    {
    }
    size_t edge_count_;
    double graph_error_;
    LoopErrorMap loop_errors_;
  };

public:
    RRRG2OWrapper(g2o::SparseOptimizer * optimizer):optimizer_(optimizer){}
    Result optimize(const EdgeIdSet& active_loops, const int n_iterations);

    // void pushOdomEdge(EdgeId & e);
    // void pushLoopEdge(EdgeId & e);
protected:
    g2o::SparseOptimizer * optimizer_;
    EdgePtrSet odom_edges_;
    EdgeIdToEdgePtrMap loop_edges_;

    // walk over all edges in graph and sort them into odometry edges and loop edges
    void parseGraph();
};

} // end of namespaces slamuk
#endif
