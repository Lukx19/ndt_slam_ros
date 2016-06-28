#ifndef GRAPH_SLAM_UK_NDT_GRID2D_HOLDER
#define GRAPH_SLAM_UK_NDT_GRID2D_HOLDER

#include <graph_slam_uk/ndt/ndt_grid2d.h>

namespace slamuk
{
class NDTGrid2DHolder : public DataHolder<NDTGrid2D>
{
public:
  explicit NDTGrid2DHolder(NDTGrid2D* data) : data_(data)
  {
  }
  Eigen::Vector2d getCentroid() override
  {
    return data_->getCentroid();
  }
  NDTGrid2D* getDataPointer() override
  {
    return data_;
  }

private:
  NDTGrid2D* data_;
}
}

#endif