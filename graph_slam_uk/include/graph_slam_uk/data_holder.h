#ifndef GRAPH_SLAM_UK_DATA_HOLDER
#define GRAPH_SLAM_UK_DATA_HOLDER

#include <Eigen/Dense>

namespace slamuk
{
template <typename Data>
class DataHolder
{
public:
  virtual Eigen::Vector2d getCentroid() = 0;
  virtual Data* getDataPointer() = 0;
}
}
#endif
