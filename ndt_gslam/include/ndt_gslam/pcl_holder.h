#ifndef NDT_GSLAM_PCL_HOLDER
#define NDT_GSLAM_PCL_HOLDER

#include <ndt_gslam/ndt/ndt_grid2d.h>
#include <pcl/point_cloud.h>

namespace slamuk
{
template <typename PointType>
class NDTGrid2DHolder
{
  typedef pcl::PointCloud<PointType> Pcl;
  typedef typename Pcl::Ptr PclPtr;
  typedef PointType Point;

public:
  explicit NDTGrid2DHolder(const PclPtr &data) : data_(data)
  {
  }
  /**
   * @brief      Gets the centroid of NDT grid.
   *
   * @return     The centroid.
   */
  Eigen::Vector2d getCentroid() const
  {
    Eigen::Vector4f centroid;
    pcl::compute3DCentroid(*data_, centroid);
    return centroid.head(2).cast<double>();
  }

  /**
   * @brief      Radius from the centroid to the farthest cell
   *
   * @return     The radius.
   */
  double getRadius() const
  {
    std::vector<float> ranges_horizontal;
    std::vector<float> ranges_vertical;

    for (auto &pt : data_.points) {
      ranges_vertical.push_back(std::abs(pt.x));
      ranges_horizontal.push_back(std::abs(pt.y));
    }
    // sort lengths 0,1,....n
    std::sort(ranges_horizontal.begin(), ranges_horizontal.end());
    std::sort(ranges_vertical.begin(), ranges_vertical.end());
    // reject 1/3 of the most distante points
    size_t idy = (ranges_vertical.size() / 3) * 2;
    size_t idx = (ranges_horizontal.size() / 3) * 2;
    return std::max(std::abs(ranges_vertical[idy]),
                    std::abs(ranges_horizontal[idx]));
  }

  /**
   * @brief      Gets the NDT grid pointer
   *
   * @return     The data.
   */
  PclPtr &getData() const
  {
    return data_;
  }

  /**
   * @brief      updates oring of the NDT grid
   *
   * @param[in]  new_pose  The new pose
   */
  void updatePosition(const Eigen::Vector3d &new_pose)
  {
    data_->setOrigin(new_pose);
  }

private:
  PclPtr data_;
};
}

#endif
