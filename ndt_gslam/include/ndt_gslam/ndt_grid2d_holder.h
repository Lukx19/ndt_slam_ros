#ifndef NDT_GSLAM_NDT_GRID2D_HOLDER
#define NDT_GSLAM_NDT_GRID2D_HOLDER

#include <ndt_gslam/ndt/ndt_grid2d.h>

namespace slamuk
{
template <typename CellType, typename PointType>
class NDTGrid2DHolder
{
public:
  typedef NDTGrid2D<CellType, PointType> Grid;
  typedef typename Grid::Ptr GridPtr;
  typedef typename Grid::ConstPtr GridConstPtr;
  typedef CellType Cell;
  typedef PointType Point;

public:
  explicit NDTGrid2DHolder(const GridPtr &data) : data_(data)
  {
  }
  /**
   * @brief      Gets the centroid of the ndt grid
   *
   * @return     The centroid.
   */
  Eigen::Vector2d getCentroid() const
  {
    return data_->getCentroid();
  }
  /**
   * @brief      Gets the radius measured from centroid
   *
   * @return     The radius.
   */
  double getRadius() const
  {
    return data_->getRadius();
  }
  const GridPtr &getData() const
  {
    return data_;
  }

  /**
   * @brief      Updates origin of the NDT grid
   *
   * @param[in]  new_pose  The new pose
   */
  void updatePosition(const Eigen::Vector3d &new_pose)
  {
    data_->setOrigin(new_pose);
  }

private:
  GridPtr data_;
};
}

#endif
