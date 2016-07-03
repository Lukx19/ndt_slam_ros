#ifndef GRAPH_SLAM_UK_NDT_MAPPER
#define GRAPH_SLAM_UK_NDT_MAPPER

#include <ros/ros.h>
#include <map>
#include <memory>
#include <Eigen/Dense>
#include <opencv/cv.h>
#include <graph_slam_uk/utils/eigen_tools.h>
#include <graph_slam_uk/ndt/ndt_grid2d.h>
#include <graph_slam_uk/ndt/ndt_grid2d_interface.h>
#include <nav_msgs/OccupancyGrid.h>
#include <graph_slam_uk/utils/msgs_conversions.h>
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
public:
  typedef NDTGrid2D<CellType, PointType> Grid;

private:
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

  const cv::Mat &getOccupancyMap() const;
  nav_msgs::OccupancyGrid calcOccupancyGridMsg() const;

private:
  size_t width_, height_;
  float resolution_;  // [m/cell]
  ros::Time map_recalc_time_;
  cv::Mat map_;
  FrameStorage grids_;

  // return id of frame if found. Otherwise exception
  typename FrameStorage::iterator getFrameIterator(const NDTGrid2DPtr &frame);
  void addToMap(const Grid &frame);
  cv::Mat frameToMat(const OccupancyGrid &occ_grid) const;
  std::pair<size_t, size_t> calcCoordinates(double x, double y, size_t width,
                                            size_t height,
                                            double cell_size) const;
};

// ////////IMPLEMENTATION////

// template <typename CellType, typename PointType>
// typename NDTMapper<CellType, PointType>::Grid *
// NDTMapper<CellType, PointType>::addFrame(const Grid &frame)
// {
//   grids_.push_back(frame);
//   addToMap(grids_.back());
//   return grids_.back().get();
// }
template <typename CellType, typename PointType>
NDTMapper<CellType, PointType>::NDTMapper()
  : width_(200), height_(200), resolution_(0.25f), map_(200, 200, CV_8UC1)
{
}

template <typename CellType, typename PointType>
void NDTMapper<CellType, PointType>::addFrame(const NDTGrid2DPtr &frame,
                                              const ros::Time &capture_time)
{
  grids_.push_back(frame);
  addToMap(*grids_.back());
  grids_.back()->setTimestamp(capture_time.toSec());
}

template <typename CellType, typename PointType>
void NDTMapper<CellType, PointType>::addFrame(NDTGrid2DPtr &&frame,
                                              const ros::Time &capture_time)
{
  grids_.push_back(std::move(frame));
  addToMap(*grids_.back());
  grids_.back()->setTimestamp(capture_time.toSec());
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
  map_ = cv::Scalar(0);
  std::sort(grids_.begin(), grids_.end(),
            [](const NDTGrid2DPtr &a,
               const NDTGrid2DPtr &b) { return a->operator<(*b); });
  for (auto &&grid : grids_) {
    addToMap(*grid);
  }
  map_recalc_time_ = calc_time;
}

template <typename CellType, typename PointType>
void NDTMapper<CellType, PointType>::addToMap(const Grid &frame)
{
  OccupancyGrid occ_grid = frame.createOccupancyGrid();
  resolution_ = occ_grid.resolution_;
  cv::Mat src = frameToMat(occ_grid);
  // cv::Mat trans_src(frame_mat.size(),frame_mat.type());
  float angle = eigt::getAngleDiffrence(Pose(0, 0, 0), frame.getOrigin());
  // get rotation matrix for rotating the image around its center
  cv::Point2f center(src.cols / 2.0f, src.rows / 2.0f);
  cv::Mat rot = cv::getRotationMatrix2D(center, angle, 1.0);
  // determine bounding rectangle
  cv::Rect bbox = cv::RotatedRect(center, src.size(), angle).boundingRect();
  // adjust transformation matrix
  rot.at<double>(0, 2) += bbox.width / 2.0 - center.x;
  rot.at<double>(1, 2) += bbox.height / 2.0 - center.y;
  cv::Mat rotated;
  cv::warpAffine(src, rotated, rot, bbox.size());
  // blend rotated frame with global map
  auto coords_pair = calcCoordinates(occ_grid.origin_(0), occ_grid.origin_(1),
                                     width_, height_, occ_grid.resolution_);
  cv::Rect roi(coords_pair.first, coords_pair.second, rotated.cols,
               rotated.rows);
  addWeighted(map_(roi), 0.3, rotated, 0.7, 0.0, map_(roi));
}

template <typename CellType, typename PointType>
nav_msgs::OccupancyGrid
NDTMapper<CellType, PointType>::calcOccupancyGridMsg() const
{
  nav_msgs::OccupancyGrid occ_msg;
  occ_msg.info.height = height_;
  occ_msg.info.width = width_;
  occ_msg.info.resolution = resolution_;
  occ_msg.info.origin = EigenToPoseMsg(Eigen::Vector3d::Zero());
  occ_msg.info.map_load_time = map_recalc_time_;
  for (int row = 0; row < map_.rows; ++row) {
    for (int col = 0; col < map_.cols; ++col) {
      if (map_.at<int>(row, col) == 0 || map_.at<int>(row, col) > 101)
        occ_msg.data.push_back(-1);
      else
        occ_msg.data.push_back(
            static_cast<uint8_t>(map_.at<int>(row, col) - 1));
    }
  }
  return occ_msg;
}

//// PROTECTED ////////
template <typename CellType, typename PointType>
cv::Mat
NDTMapper<CellType, PointType>::frameToMat(const OccupancyGrid &occ_grid) const
{
  cv::Mat frame_mat(occ_grid.height_, occ_grid.width_, CV_8UC1);
  for (size_t i = 0; i < occ_grid.cells_.size(); ++i) {
    if (occ_grid.cells_[i] == -1 || occ_grid.cells_[i] > 100)
      frame_mat.data[i] = 0;
    else
      frame_mat.data[i] = occ_grid.cells_[i] + 1;
  }
}

template <typename CellType, typename PointType>
typename NDTMapper<CellType, PointType>::FrameStorage::iterator
NDTMapper<CellType, PointType>::getFrameIterator(const NDTGrid2DPtr &frame)
{
  auto iter = grids_.find(frame);
  if (iter == grids_.end())
    throw std::invalid_argument("Id of frame doesn't exist");
  return iter;
}

template <typename CellType, typename PointType>
std::pair<size_t, size_t> NDTMapper<CellType, PointType>::calcCoordinates(
    double x, double y, size_t width, size_t height, double cell_size) const
{
  std::pair<size_t, size_t> coords;
  double minx = cell_size * (width / 2);
  double maxy = cell_size * (height / 2);
  coords.first =
      static_cast<size_t>(std::floor((x - minx + cell_size / 2) / cell_size));
  coords.second =
      static_cast<size_t>(std::floor((-y + maxy + cell_size / 2) / cell_size));
  return coords;
}

}  // end of namespace slamuk
#endif
