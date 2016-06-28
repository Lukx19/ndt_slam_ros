#include <graph_slam_uk/rrr_clusterizer.h>

namespace slamuk
{
namespace internal
{
std::vector<size_t> Clusterizer::clusterizeBulk(const IntPairSet & loops)
{
  std::vector<size_t> new_clusters;
  closed_clusters_.clear();
  updateClusters(loops);
  if(closed_clusters_.size() != 0)
    std::cerr<<"Something wrong: Update should not have filled anything to closed loops\n";
  for(size_t i  =0; i< active_clusters_.size(); ++i){
    for(auto && loop:active_clusters_[i].loops_){
      loop_to_id_map_.insert(std::make_pair(loop,i));
    }

    closed_clusters_.push_back(std::move(active_clusters_[i]));
    new_clusters.emplace_back(i);
  }
  return new_clusters;
}

// return vector of new closed clusters ids
std::vector<size_t> Clusterizer::updateClusters(const IntPairSet & loops)
{
  std::vector<size_t> recent_closed;
  for(auto && edge:loops){
    clusterizeOne(edge);
  }
  // decrease time to live for all current clusters
  for(auto && cluster: active_clusters_){
    if(cluster.ttl_ == 0){
      // add move semantics in future
      closed_clusters_.push_back(cluster);
      size_t id  = closed_clusters_.size()-1;
      recent_closed.push_back(id);
      // adding to map of loops
      for(auto && loop:cluster.loops_){
        loop_to_id_map_.insert(std::make_pair(loop,id));
      }
    }
    --cluster.ttl_;
  }
  // remove old closed clusters
  active_clusters_.erase(std::remove_if(active_clusters_.begin(), active_clusters_.end(),
                 [](const Cluster & edge) {
                  if(edge.ttl_ < 0)
                    return true;
                  return false; }),active_clusters_.end());
  return recent_closed;
}


void Clusterizer::clusterizeOne(const IntPair& loop)
{
  int start = std::max(loop.first, loop.second);
  int end = std::min(loop.first, loop.second);
  bool found_cluster =  false;
  if (active_clusters_.empty()) {
    active_clusters_.emplace_back(Cluster(loop));
  } else {
    // Search for a cluster where it can belong
    for (size_t i = 0; i < active_clusters_.size(); i++) {
      if (active_clusters_[i].contains(start, end, threshold_)) {
        active_clusters_[i].updateClusterVals(loop);
        found_cluster = true;
        break;
      }
    }
    if (!found_cluster)  // Does not belong to any existing cluster.
    {
      active_clusters_.emplace_back(Cluster(loop));
    }
  }

#if 0
        if(0)
        {
            std::cout<<" \% Clusters formed "<<_clustersFound.size()<<std::endl;
            std::cout<<"limits = [ "<<std::endl;
            for(size_t i=0 ; i< _clustersFound.size() ; i++)
            {
                std::cout<<i<<" -> sz "<<_clustersFound[i].size<<" :: ";
                std::cout<<" "<<_clustersFound[i].startLow<<" "<<_clustersFound[i].startHigh<<" ";
                std::cout<<" "<<_clustersFound[i].endLow<<" "<<_clustersFound[i].endHigh<<std::endl;;

            }
            std::cout<<std::endl;
            std::cout<<"]; "<<std::endl;


            std::cout<<"membership =[ ";
            for(size_t i=0; i<membership.size();i++)
            {
                std::cout<<membership[i]<<" ";
            }
            std::cout<<std::endl;
            std::cout<<"]; "<<std::endl;
        }
#endif
}

Cluster & Clusterizer::getCluster(size_t id){
    assert(id > closed_clusters_.size()-1);
    return closed_clusters_[id];
}

const Cluster& Clusterizer::getCluster(size_t id) const
{
  assert(id > closed_clusters_.size() - 1);
  return closed_clusters_[id];
}

size_t Clusterizer::getClusterID(const IntPair & loop){
  return loop_to_id_map_[loop];
}

size_t Clusterizer::clusterCount()
{
  return closed_clusters_.size();
}

bool Clusterizer::deleteCluster(size_t id)
{
  auto it = closed_clusters_.begin();
  std::advance(it, id);
  std::for_each(
      it->loops_.begin(), it->loops_.end(),
      [&](const IntPair& loop) { loop_to_id_map_.erase(loop); });
  closed_clusters_.erase(it);
  return true;
}
}} // end of namespaces
