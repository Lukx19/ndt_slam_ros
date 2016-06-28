#ifndef GRAPH_SLAM_UK_RRR_CLUSTERIZER
#define GRAPH_SLAM_UK_RRR_CLUSTERIZER

#include <vector>
#include <cmath>
#include <set>
#include <map>
#include <deque>
#include <iostream>
#include <algorithm>
#include <cassert>

namespace slamuk
{
namespace internal
{
struct Cluster
{
public:
  typedef std::pair<int,int> Edge;
  int start_low_, start_high_;
  int end_low_, end_high_;
  int size_;
  bool valid_;
  int ttl_;  // how many iterations of updates it should be in active clusters
  std::set<Edge> loops_;

  Cluster()
    : start_low_(-1)
    , start_high_(-1)
    , end_low_(-1)
    , end_high_(-1)
    , size_(0)
    , valid_(false)
    , ttl_(-1)
  {
  }
  Cluster(int start, int end)
    : start_low_(start)
    , start_high_(start)
    , end_low_(end)
    , end_high_(end)
    , size_(1)
    ,valid_(true)
    , ttl_(1)
  {
  }
  Cluster(const Edge & loop)
  {
    Cluster(loop.first, loop.second);
    loops_.insert(loop);
  }

  bool contains(int start, int end, int threshold)
  {
    return (std::abs(start - start_high_) < threshold or
            std::abs(start - start_low_) < threshold) and
           (std::abs(end - end_high_) < threshold or
            std::abs(end - end_low_) < threshold);
  }

  void updateClusterVals(int start, int end){
    if (start < start_low_)
      start_low_ = start;
    if (start > start_high_)
      start_high_ = start;

    if (end < end_low_)
      end_low_ = end;
    if (end > end_high_)
      end_high_ = end;
    ++size_;
  }
  void updateClusterVals(const Edge & edge){
    updateClusterVals(edge.first,edge.second);
    loops_.insert(edge);
  }
};

class Clusterizer
{
public:
  typedef std::pair<int, int> IntPair;
  typedef std::set<IntPair> IntPairSet;

private:
  typedef std::map<IntPair, int> LoopToClusterIDMap;
  // typedef IntPairIDMap LoopToClusterIDMap;
  // typedef IDintPairSetMap ClusterIDtoLoopsMap;
  typedef std::map<int, IntPairSet> ClusterIDtoLoopsMap;

public:
  Clusterizer():threshold_(50){}
  Clusterizer(int treshold):threshold_(treshold){}

  std::vector<size_t> clusterizeBulk(const IntPairSet & loops);
  std::vector<size_t> updateClusters(const IntPairSet & loops);
  Cluster & getCluster(size_t id);
  const Cluster & getCluster(size_t id) const;
  size_t getClusterID(const IntPair&);
  size_t clusterCount();
  bool deleteCluster(size_t id);

protected:
  int top_cluster_id_;
  int threshold_;
  std::vector<Cluster> closed_clusters_;
  std::deque<Cluster> active_clusters_;
  std::map<IntPair,size_t> loop_to_id_map_;

  void clusterizeOne(const IntPair& loop);


  std::vector<Cluster> _clustersFound; // holds information about clusters
  ClusterIDtoLoopsMap clusterIDtoLoopsMap;  // holds all clusters (cluster = set
                                            // of loop closure edges) by id
  LoopToClusterIDMap loopToClusterIDMap;    // key = loop closure val = index of
                                            // cluster
};
}  // end of internal namespace
}  // end of slamuk
#endif
