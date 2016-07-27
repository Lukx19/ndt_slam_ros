#ifndef GRAPH_SLAM_UK_NDT_MAPPER
#define GRAPH_SLAM_UK_NDT_MAPPER

#include <graph_slam_uk/ndt/ndt_grid2d.h>
#include <graph_slam_uk/utils/eigen_tools.h>
#include <graph_slam_uk/utils/msgs_conversions.h>
#include <nav_msgs/OccupancyGrid.h>
#include <pcl/point_cloud.h>
#include <ros/ros.h>
#include <Eigen/Dense>
#include <map>
#include <memory>

/*
Grid needs to have implemented:
operator<
setOrigin(Pose)
setTimestamp(double)
getTimestamp(double)

*/
namespace slamuk
{
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
  void addFrame(const NDTGrid2DPtr &frame, const ros::Time &capture_time);
  void addFrame(NDTGrid2DPtr &&frame, const ros::Time &capture_time);
  void removeFrame(const NDTGrid2DConstPtr &frame);
  void updateFrame(const NDTGrid2DPtr &frame, const Pose &new_pose);
  NDTGrid2DPtr combineFrames(const NDTGrid2DPtr &a, const NDTGrid2DPtr &b);
  void recalc(const ros::Time &calc_time);

  nav_msgs::OccupancyGrid calcOccupancyGridMsg() const;

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

//// PROTECTED ////////

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
