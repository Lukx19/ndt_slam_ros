#include <graph_slam_uk/rrr_g2o_wrapper.h>

using namespace slamuk;

//template<typename VertexType,typename EdgeType>
 RRRG2OWrapper::Result RRRG2OWrapper/*<VertexType,EdgeType>*/::optimize(const EdgeIdSet & active_loops,
                const int n_iterations)
{
    parseGraph();
    // initialize all nodes needed for optimalization
    EdgePtrSet  active_edges(odom_edges_.begin(),odom_edges_.end());
    for(auto && loop:active_loops){

        active_edges.insert(loop_edges_[loop]);
    }

    optimizer_->setVerbose(false);
    optimizer_->findGauge()->setFixed(true);
    optimizer_->initializeOptimization(active_edges);
    optimizer_->optimize(n_iterations,false);
    optimizer_->computeActiveErrors();
    optimizer_->push();
    // output calculated values
    Result res;
    res.loop_errors_.clear();
    for(auto && loop_id : active_loops){
//        EdgeType * e = dynamic_cast<EdgeType *>(loop_edges_[loop_id]);
        res.loop_errors_.emplace(std::make_pair(loop_id,loop_edges_[loop_id]->chi2()));
    }
    res.edge_count_ = optimizer_->activeEdges().size();
    res.graph_error_ = optimizer_->activeChi2();
    optimizer_->pop();
    return res;
}

//template<typename VertexType,typename EdgeType>
void RRRG2OWrapper/*<VertexType,EdgeType>*/::parseGraph()
{
    g2o::OptimizableGraph::EdgeSet::iterator
            it = optimizer_->edges().begin(),
            e_it = optimizer_->edges().end();

    for (;it != e_it; ++it) {
    // for(const g2o::HyperGraph::Edge * edge : optimizer_->edges()){
        int from = (*it)->vertex(0)->id();
        int to = (*it)->vertex(1)->id();
        if(std::abs(to-from) > 1 ){
            // loop closure
            loop_edges_.emplace(std::make_pair(std::make_pair(from,to),*it));
        }else{
            // odom
            odom_edges_.insert(*it);
        }
    }
}

