#ifndef NDT_GSLAM_NDT_MAPPER
#define NDT_GSLAM_NDT_MAPPER

#include <nav_msgs/OccupancyGrid.h>
#include <ndt_gslam/ndt/ndt_grid2d.h>
#include <ndt_gslam/utils/eigen_tools.h>
#include <ndt_gslam/utils/msgs_conversions.h>
#include <pcl/point_cloud.h>
#include <ros/ros.h>
#include <Eigen/Dense>
#include <map>
#include <memory>

namespace slamuk
{
/**
 * @brief      {The class merging individual ndt_grids into one map}
 *
 * @tparam     CellType   { type of NDT cell}
 * @tparam     PointType  { type of PCL point used for initialization of
 * ndt_grid}
 */
template <typename CellType, typename PointType>
class NDTMapper
{
private:
  typedef NDTGrid2D<CellType, PointType> Grid;
  typedef typename pcl::PointCloud<PointType> Pcl;

public:
  typedef typename Grid::Ptr NDTGrid2DPtr;
  typedef typename Grid::ConstPtr NDTGrid2DConstPtr;
  typedef typename Grid::Pose Pose;
  typedef std::vector<NDTGrid2DPtr> FrameStorage;

public:
  NDTMapper();

  /**
   * @brief      Adds a frame to the map. Map is updated.
   *
   * @param[in]  frame         The frame
   * @param[in]  capture_time  The capture time
   */
  void addFrame(const NDTGrid2DPtr &frame, const ros::Time &capture_time);

  /**
   * @brief      Adds a frame frame to the map. Map is updated.
   *
   * @param[in]  frame    The frame
   * @param[in]  capture_time  The capture time
   */
  void addFrame(NDTGrid2DPtr &&frame, const ros::Time &capture_time);

  /**
   * @brief      Removes a frame.
   *
   * @param[in]  frame  The frame
   */
  void removeFrame(const NDTGrid2DConstPtr &frame);

  /**
   * @brief      Changes frames origin. Do not recalculate the whole map.
   *
   * @param[in]  frame     The frame
   * @param[in]  new_pose  The new pose
   */
  void updateFrame(const NDTGrid2DPtr &frame, const Pose &new_pose);

  /**
   * @brief      Merge two frames into one.
   *
   * @param[in]  a     The first frame
   * @param[in]  b     The second frame
   *
   * @return     The pointer to new created frame used in map.
   */
  NDTGrid2DPtr combineFrames(const NDTGrid2DPtr &a, const NDTGrid2DPtr &b);

  /**
   * @brief      Redraws whole map.
   *
   * @param[in]  calc_time  The calculate time of recalculation
   */
  void recalc(const ros::Time &calc_time);

  /**
   * @brief      Calculates the occupancy grid message.
   *
   * @return     The occupancy grid message.
   */
  nav_msgs::OccupancyGrid calcOccupancyGridMsg() const;

  /**
   * @brief      Gets the PCL map representation.
   *
   * @return     The PCL map.
   */
  typename Pcl::Ptr getPclMap() const
  {
    return means_;
  }

private:
  size_t width_, height_;
  float resolution_;  // [m/cell]
  ros::Time map_recalc_time_;
  OccupancyGrid map_;
  Grid map_ndt_;
  FrameStorage grids_;
  typename Pcl::Ptr means_;
  // return id of frame if found. Otherwise exception
  typename FrameStorage::iterator getFrameIterator(const NDTGrid2DPtr &frame);
  void addToMap(const Grid &frame);
};

// ////////IMPLEMENTATION////

template <typename CellType, typename PointType>
NDTMapper<CellType, PointType>::NDTMapper()
  : width_(200)
  , height_(200)
  , resolution_(0.25f)
  , map_()
  , map_ndt_(Eigen::Vector3d(0, 0, 0))
  , grids_()
  , means_(new Pcl())
{
  map_ndt_.setCellSize(resolution_);
}

template <typename CellType, typename PointType>
void NDTMapper<CellType, PointType>::addFrame(const NDTGrid2DPtr &frame,
                                              const ros::Time &capture_time)
{
  grids_.push_back(frame);
  addToMap(*grids_.back());
  grids_.back()->setTimestamp(capture_time.toSec());
  map_ = map_ndt_.createOccupancyGrid();
  means_ = map_ndt_.getMeansTransformed();
  std::cout << "node added" << std::endl;
}

template <typename CellType, typename PointType>
void NDTMapper<CellType, PointType>::addFrame(NDTGrid2DPtr &&frame,
                                              const ros::Time &capture_time)
{
  grids_.push_back(std::move(frame));
  addToMap(*grids_.back());
  grids_.back()->setTimestamp(capture_time.toSec());
  map_ = map_ndt_.createOccupancyGrid();
  means_ = map_ndt_.getMeansTransformed();
  std::cout << "node added" << std::endl;
}

template <typename CellType, typename PointType>
void NDTMapper<CellType, PointType>::removeFrame(const NDTGrid2DConstPtr &frame)
{
  auto frame_iter = getFrameIterator(frame);
  grids_.erase(frame_iter);
}

template <typename CellType, typename PointType>
void NDTMapper<CellType, PointType>::updateFrame(const NDTGrid2DPtr &frame,
                                                 const Pose &pose)
{
  getFrameIterator(frame)->setOrigin(pose);
}

template <typename CellType, typename PointType>
typename NDTMapper<CellType, PointType>::NDTGrid2DPtr
NDTMapper<CellType, PointType>::combineFrames(const NDTGrid2DPtr &a,
                                              const NDTGrid2DPtr &b)
{
  Grid *frame_a = getFrameIterator(a)->get();
  Grid *frame_b = getFrameIterator(b)->get();
  frame_a->mergeInTraced(*frame_b, true, true);
  frame_a->setTimestamp(frame_b->getTimestamp());
  removeFrame(frame_b);
}

template <typename CellType, typename PointType>
void NDTMapper<CellType, PointType>::recalc(const ros::Time &calc_time)
{
  std::sort(grids_.begin(), grids_.end(),
            [](const NDTGrid2DPtr &a, const NDTGrid2DPtr &b) {
              return a->operator<(*b);
            });
  map_ndt_.clear();
  for (const NDTGrid2DPtr &grid : grids_) {
    addToMap(*grid);
  }
  map_ = map_ndt_.createOccupancyGrid();
  means_ = map_ndt_.getMeansTransformed();
  map_recalc_time_ = calc_time;
}

template <typename CellType, typename PointType>
void NDTMapper<CellType, PointType>::addToMap(const Grid &frame)
{
  map_ndt_.mergeInTraced(frame, true, true);
}

template <typename CellType, typename PointType>
nav_msgs::OccupancyGrid
NDTMapper<CellType, PointType>::calcOccupancyGridMsg() const
{
  nav_msgs::OccupancyGrid occ_msg;
  occ_msg.info.height = map_.height_;
  occ_msg.info.width = map_.width_;
  occ_msg.info.resolution = map_.resolution_;
  occ_msg.info.origin = EigenToPoseMsg(map_.origin_);
  occ_msg.info.map_load_time = map_recalc_time_;
  occ_msg.data = map_.cells_;
  return occ_msg;
}

//// PRIVATE ////////

template <typename CellType, typename PointType>
typename NDTMapper<CellType, PointType>::FrameStorage::iterator
NDTMapper<CellType, PointType>::getFrameIterator(const NDTGrid2DPtr &frame)
{
  auto iter = grids_.find(frame);
  if (iter == grids_.end())
    throw std::invalid_argument("Id of frame doesn't exist");
  return iter;
}

}  // end of namespace slamuk
#endif
