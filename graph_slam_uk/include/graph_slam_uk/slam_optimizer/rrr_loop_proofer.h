// RRR - Robust Loop Closing over Time
// Copyright (C) 2014 Y.Latif, C.Cadena, J.Neira
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef GRAPH_SLAM_UK_RRR_LOOP_PROOFER
#define GRAPH_SLAM_UK_RRR_LOOP_PROOFER

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <boost/math/distributions/chi_squared.hpp>
#include <graph_slam_uk/slam_optimizer/rrr_clusterizer.h>

namespace slamuk
{
template <class Optimizer>
class RRRLoopProofer
{
private:
  typedef std::vector<size_t> IdVector;
  typedef std::pair<int, int> VerticesIntPair;
  typedef std::set<VerticesIntPair> VerticesIntPairSet;

public:
  typedef std::pair<size_t, size_t> VerticesPair;
  typedef std::set<VerticesPair> VerticesPairSet;

public:
  RRRLoopProofer(Optimizer* opt, int clustering_treshold = 50,
                 int numIterations = 4, int dimmention = 2)
    : optimizer_(opt)
    , clustering_treshold_(clustering_treshold)
    , n_iterations_(numIterations)
    , edge_dimention_(dimmention)
  {
  }

  bool optimizeInc(const VerticesPairSet& possible_loops,
                   VerticesPairSet& bad_loops, VerticesPairSet& good_loops);

protected:
  // TODO: remove when done
  Optimizer* optimizer_;
  int clustering_treshold_;
  int n_iterations_;
  int edge_dimention_;

  // used fo incremental adding
  IdVector last_good_clusters_;
  internal::Clusterizer clusterizer_;

  bool intraClusterConsistent(size_t cluster_id, VerticesIntPairSet& bad_loops);
  VerticesIntPairSet gatherEdges(const IdVector& clusterList);
  bool interClusterConsistent(IdVector& H, IdVector& goodSet,
                              IdVector& RejectSet);
  void optimizeCluster(size_t cluster_id, VerticesIntPairSet& bad_loops,
                       VerticesIntPairSet& good_loops);
  double chi2(int dof)
  {
    if (dof > 0)
      return boost::math::quantile(boost::math::chi_squared(dof), 0.95);
    else {
      std::cerr << __LINE__ << " dof <= 0 ? " << std::endl;
      return 0;
    }
  }
};

//////////////////////////////////////IMPLEMENTATION //////////////////////
template <class Optimizer>
bool RRRLoopProofer<Optimizer>::intraClusterConsistent(
    size_t cluster_id, VerticesIntPairSet& bad_loops)
{
  typename Optimizer::Result res;
  res = optimizer_->optimize(clusterizer_.getCluster(cluster_id).loops_,
                             n_iterations_);

  if (res.graph_error_ < chi2(edge_dimention_ * res.edge_count_)) {
    // IntPairDoubleMap::iterator eIt = chi2LinkErrors.begin(),
    //                            eEnd = chi2LinkErrors.end();
    // VerticesIntPairSet::iterator it = currentCluster.begin(), end =
    // currentCluster.end();
    for (auto&& edge_err : res.loop_errors_) {
      if (edge_err.second > chi2(edge_dimention_)) {
        bad_loops.insert(edge_err.first);
        clusterizer_.getCluster(cluster_id).loops_.erase(edge_err.first);
      }
    }
    // for (; eIt != eEnd; eIt++) {
    //   if (eIt->second > chi2(edge_dimention_)) {
    //     clusterizer_.setClusterID(eIt->first, ID_IGNORE);
    //   }
    if (!clusterizer_.getCluster(cluster_id).loops_.empty()) {
      std::cerr << " Cluster " << cluster_id << " survived with"
                << clusterizer_.getCluster(cluster_id).size_ << " links "
                << std::endl;
      return true;
    } else {
      clusterizer_.getCluster(cluster_id).valid_ = false;
      std::cerr << " Cluster " << cluster_id << " has been emptied!"
                << std::endl;
      return false;
    }

  } else {
    clusterizer_.getCluster(cluster_id).valid_ = false;
    clusterizer_.getCluster(cluster_id).loops_.clear();
    std::cerr << " Cluster " << cluster_id << " has been eliminated!"
              << std::endl;
    bad_loops.insert(clusterizer_.getCluster(cluster_id).loops_.begin(),
                     clusterizer_.getCluster(cluster_id).loops_.end());
    return false;
  }

  return false;
}

template <class Optimizer>
typename RRRLoopProofer<Optimizer>::VerticesIntPairSet
RRRLoopProofer<Optimizer>::gatherEdges(const IdVector& cluster_ids)
{
  VerticesIntPairSet links;
  for (size_t id : cluster_ids) {
    links.insert(clusterizer_.getCluster(id).loops_.begin(),
                 clusterizer_.getCluster(id).loops_.end());
  }
  return links;
}

template <class Optimizer>
bool RRRLoopProofer<Optimizer>::interClusterConsistent(IdVector& candidates,
                                                       IdVector& good_clusters,
                                                       IdVector& bad_clusters)
{
  if (candidates.empty())
    return true;

  VerticesIntPairSet activeLinks(gatherEdges(good_clusters));
  auto temp = gatherEdges(candidates);
  std::move(temp.begin(), temp.end(),
            std::inserter(activeLinks, activeLinks.begin()));

  typename Optimizer::Result res;
  res = optimizer_->optimize(activeLinks, n_iterations_);

  double activeChi2Graph = res.graph_error_;
  size_t activeEdgeCount = res.edge_count_;

  double allLinksError = 0;
  for (auto&& edge_err : res.loop_errors_) {
    allLinksError += edge_err.second;
  }

  if (activeChi2Graph < chi2(edge_dimention_ * activeEdgeCount) and
      allLinksError < chi2(edge_dimention_ * activeLinks.size())) {
    std::copy(candidates.begin(), candidates.end(),
              std::back_inserter(good_clusters));
    return true;  // all done .. we found a consistent solution
  } else {
    // Find which cluster is causing the problems
    // Iterate over everything
    // Find the error for each cluster
    // Sum cluster-wise
    std::map<size_t, double> error_map;
    for (const auto& edge_err : res.loop_errors_) {
      error_map[clusterizer_.getClusterID(edge_err.first)] += edge_err.second;
    }

    double min_CI = 0;
    long rejectID = -1;
    for (size_t cluster_id : candidates) {
      double CI = error_map[cluster_id] /
                  chi2(edge_dimention_ *
                       clusterizer_.getCluster(cluster_id).loops_.size());
      if (CI >= min_CI)  // Just looking for the ones in candidates
      {
        min_CI = CI;
        rejectID = static_cast<long>(cluster_id);
      }
    }
    candidates.erase(candidates.begin() + rejectID);
    bad_clusters.push_back(static_cast<size_t>(rejectID));

    return interClusterConsistent(candidates, good_clusters, bad_clusters);
  }
  return true;
}

// returns true if some optimalization was done
template <class Optimizer>
bool
RRRLoopProofer<Optimizer>::optimizeInc(const VerticesPairSet& possible_loops,
                                       VerticesPairSet& bad_loops_fin,
                                       VerticesPairSet& good_loops_fin)
{
  bad_loops_fin.clear();
  good_loops_fin.clear();
  // compatibility change from pair of size_t to pair of int
  VerticesIntPairSet possible_loops_int;
  for (auto&& loop : possible_loops) {
    possible_loops_int.insert(std::make_pair(static_cast<int>(loop.first),
                                             static_cast<int>(loop.second)));
  }
  // compatibility of output
  VerticesIntPairSet bad_loops;
  VerticesIntPairSet good_loops;

  // add possible loops to clusters
  std::vector<size_t> new_clusters =
      clusterizer_.updateClusters(possible_loops_int);
  std::cout << "Added to clusters:\n";
  for (auto&& cluster_id : new_clusters) {
    std::cout << cluster_id;
  }
  std::cout << std::endl;
  // std::vector<size_t> valid_clusters;
  // if(new_clusters.empty()){
  //   // nothing to optimize
  //   std::cout<<"No loop closure clusters to optimize\n";
  //   return false;
  // }
  // // check for internal consistency of every cluster
  // for(auto && cluster_id : new_clusters){
  //   if(intraClusterConsistent(cluster_id,bad_loops)){
  //     valid_clusters.push_back(cluster_id);
  //   }else{
  //     clusterizer_.deleteCluster(cluster_id);
  //   }
  // }
  // // valid_clusters hold only clusters with good edges
  // // all rejected edges are added to bad_loops

  // IdVector rejected_clusters;
  // IdVector candidates;
  // std::sort(valid_clusters.begin(),valid_clusters.end());
  // while(true){
  //   IdVector current_good_clusters;
  //   std::sort(rejected_clusters.begin(),rejected_clusters.end());
  //   std::set_difference(valid_clusters.begin(), valid_clusters.end(),
  //                       rejected_clusters.begin(), rejected_clusters.end(),
  //                       std::back_inserter(current_good_clusters));
  //   VerticesIntPairSet current_loops = gatherEdges(current_good_clusters);
  //   typename Optimizer::Result res;
  //   res = optimizer_->optimize(current_loops,n_iterations_);
  //   candidates.clear();
  //   for(auto && edge_err : res.loop_errors_){
  //     size_t id = clusterizer_.getClusterID(edge_err.first);
  //     if(!clusterizer_.getCluster(id).valid_)
  //       continue;
  //     if (edge_err.second < chi2(edge_dimention_)) {
  //       // check if not already included in candidates
  //       if(std::find(candidates.begin(),candidates.end(),id) ==
  //       candidates.end())
  //         candidates.push_back(id);
  //     }
  //   }
  //   // no more clusters to examine
  //   if(candidates.empty())
  //     break;
  //   // candidates have all new good clusters to test
  //   // merge last clusters with new clusters
  //   IdVector all_candidates;
  //   std::sort(candidates.begin(),candidates.end());
  //   std::sort(last_good_clusters_.begin(),last_good_clusters_.end());
  //   std::set_union(candidates.begin(), candidates.end(),
  //                  last_good_clusters_.begin(), last_good_clusters_.end(),
  //                  std::back_inserter(all_candidates));

  //   size_t good_clusters_size = last_good_clusters_.size();

  //   IdVector temp_rejected;
  //   interClusterConsistent(all_candidates, last_good_clusters_,
  //   temp_rejected);
  //   candidates.clear();
  //   if (last_good_clusters_.size() > good_clusters_size) {
  //          rejected_clusters.clear();
  //   }
  //   rejected_clusters.insert(rejected_clusters.end(), temp_rejected.begin(),
  //                            temp_rejected.end());
  // }
  // // all loops from rejected clusters should be outputed as rejected
  // for (size_t cluster_id : rejected_clusters) {
  //   internal::Cluster& cluster = clusterizer_.getCluster(cluster_id);
  //   cluster.valid_ = false;
  //   std::move(cluster.loops_.begin(), cluster.loops_.end(),
  //             std::inserter(bad_loops, bad_loops.begin()));
  //   clusterizer_.deleteCluster(cluster_id);
  // }
  // // all good loops in last_good_cluster should be outputed as good one
  // for (size_t cluster_id : last_good_clusters_) {
  //   internal::Cluster& cluster = clusterizer_.getCluster(cluster_id);
  //   std::copy(cluster.loops_.begin(), cluster.loops_.end(),
  //            std::inserter(good_loops, good_loops.begin()));
  // }

  // output edges in correct form
  for (auto&& loop : bad_loops) {
    bad_loops_fin.insert(std::make_pair(static_cast<size_t>(loop.first),
                                        static_cast<size_t>(loop.second)));
  }

  for (auto&& loop : good_loops) {
    good_loops_fin.insert(std::make_pair(static_cast<size_t>(loop.first),
                                         static_cast<size_t>(loop.second)));
  }
  return true;
}

// /** bad_loops[out], good_loops[in/out]*/
// template <class Optimizer>
// void RRRLoopProofer<Optimizer>::optimizeCluster(size_t cluster_id,
//                                                 VerticesIntPairSet&
//                                                 bad_loops,
//                                                 VerticesIntPairSet&
//                                                 good_loops)
// {

// }
// template <class Optimizer>
// bool RRRLoopProofer<Optimizer>::robustify(bool eraseIncorrectLinks)
// {
//   if (!gWrapper->isInitialized()) {
//     std::cerr << " Please read in a graph file with read() or set the "
//                  "optimizer with setOptimizer()" << std::endl;
//     return false;
//   }
//   // First look for intra-cluster consistency

//   IntSet consistentClusters, hypotheses, goodSet, rejectSet, tempRejectSet;

//   std::cout << "Number of Clusters found : " << clusterizer.clusterCount()
//             << std::endl;
//   std::cout << "Checking Intra cluster consistency : " << std::endl;
//   for (size_t i = 0; i < clusterizer.clusterCount(); i++) {
//     std::cout << i << " ";
//     std::cout.flush();
//     if (intraClusterConsistent(i)) {
//       consistentClusters.insert(i);
//     }
//   }
//   std::cout << "done" << std::endl << std::endl;

//   /// ------------- consistentCluster are self-consistent --------- ////

//   bool done = false;
//   while (!done) {
//     done = true;

//     IntSet::iterator cIt = consistentClusters.begin(),
//                      cEnd = consistentClusters.end();

//     for (; cIt != cEnd; cIt++) {
//       if (goodSet.find(*cIt) == goodSet.end() and
//           rejectSet.find(*cIt) == rejectSet.end())  // We can't ignore this
//                                                     // because it is nether
//                                                     in
//                                                     // goodSet nor in
//                                                     // rejectSet at the
//                                                     moment
//       {
//         hypotheses.insert(*cIt);
//       }
//     }

//     VerticesIntPairSet activeLoops;
//     gatherLinks(hypotheses, activeLoops);

//     IntPairDoubleMap linkErrors =
//         gWrapper->optimize(activeLoops, nIterations);

//     hypotheses.clear();
//     for (IntPairDoubleMap::iterator it = linkErrors.begin(),
//                                     end = linkErrors.end();
//          it != end; it++) {
//       if (it->first.first < 0)
//         continue;
//       if (it->second < chi2(edge_dimention_)) {
//         hypotheses.insert(clusterizer.getClusterID(it->first));
//         done = false;
//       }
//     }

//     size_t goodSetSize = goodSet.size();

//     interClusterConsistent(hypotheses, goodSet, tempRejectSet);

//     hypotheses.clear();

//     if (goodSet.size() > goodSetSize) {
//       rejectSet.clear();
//     }
//     rejectSet.insert(tempRejectSet.begin(), tempRejectSet.end());
//     tempRejectSet.clear();
//   }
//   std::cout << " GoodSet :";
//   for (IntSet::iterator it = goodSet.begin(), end = goodSet.end(); it != end;
//        it++) {
//     std::cout << (*it) << " ";
//   }
//   std::cout << std::endl;

//   _goodSet.insert(goodSet.begin(), goodSet.end());

//   if (eraseIncorrectLinks)
//     removeIncorrectLoops();

//   return true;
// }

// template <class Optimizer>
// bool RRRLoopProofer<Optimizer>::write(const char* filename)
// {
//   VerticesIntPairSet correctLinks;
//   gatherLinks(_goodSet, correctLinks);
//   gWrapper->write(filename, correctLinks);

//   return true;
// }

// template <class Optimizer>
// bool RRRLoopProofer<Optimizer>::removeIncorrectLoops()
// {
//   size_t count = clusterizer.clusterCount();

//   IntSet rejected;

//   for (size_t i = 0; i < count; i++) {
//     if (_goodSet.find(i) == _goodSet.end()) {
//       rejected.insert(i);
//     }
//   }
//   rejected.insert(ID_IGNORE);

//   VerticesIntPairSet rejectedLinks;

//   gatherLinks(rejected, rejectedLinks);
//   gWrapper->removeEdges(rejectedLinks);

//   // Clear up the clusters from Clusterizer

//   for (IntSet::const_iterator it = rejected.begin(), end = rejected.end();
//        it != end; it++) {
//     clusterizer.deleteCluster(*it);
//   }

//   return true;
// }

}  // end of slamuk namespace
#endif
