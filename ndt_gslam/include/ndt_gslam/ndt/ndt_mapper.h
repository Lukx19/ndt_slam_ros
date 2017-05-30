#ifndef NDT_GSLAM_NDT_MAPPER
#define NDT_GSLAM_NDT_MAPPER

#include <nav_msgs/OccupancyGrid.h>
#include <ndt_gslam/ndt/ndt_grid2d.h>
#include <ndt_gslam/utils/eigen_tools.h>
#include <ndt_gslam/utils/msgs_conversions.h>
#include <pcl/common/time.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <ros/ros.h>
#include <Eigen/Dense>
#include <map>
#include <memory>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

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
  NDTMapper(float resulution);

  /**
   * @brief      Adds a frame to the map. Map is updated.
   *
   * @param[in]  frame         The frame
   * @param[in]  capture_time  The capture time
   */
  void addFrame(const NDTGrid2DPtr &frame, const ros::Time &capture_time);

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

  NDTGrid2DConstPtr getNDTMap() const
  {
    return map_ndt_;
  }

private:
  const static constexpr float RESOLUTION_DIVISOR = 5;
  size_t width_, height_;
  float resolution_;  // [m/cell]
  ros::Time map_recalc_time_;
  OccupancyGrid map_;
  NDTGrid2DPtr map_ndt_;
  FrameStorage grids_;
  typename Pcl::Ptr means_;

  // return id of frame if found. Otherwise exception
  typename FrameStorage::iterator getFrameIterator(const NDTGrid2DPtr &frame);
  void addToMap(const Grid &frame);
  void repareOccupancyGrid();

  cv::Mat frameToMat(const OccupancyGrid &occ_grid) const;
  OccupancyGrid matToFrame(const cv::Mat &img,
                           const OccupancyGrid &occ_grid) const;
};

// ////////IMPLEMENTATION////

template <typename CellType, typename PointType>
NDTMapper<CellType, PointType>::NDTMapper(float resulution)
  : width_(200)
  , height_(200)
  , resolution_(resulution)
  , map_()
  , map_ndt_(new Grid(resolution_, Eigen::Vector3d(0, 0, 0)))
  , grids_()
  , means_(new Pcl())
{
}

template <typename CellType, typename PointType>
void NDTMapper<CellType, PointType>::addFrame(const NDTGrid2DPtr &frame,
                                              const ros::Time &capture_time)
{
  pcl::io::savePCDFile("frame" + std::to_string(grids_.size()) + ".pcd",
                       *(frame->getMeans()));
  grids_.push_back(frame);
  addToMap(*grids_.back());
  grids_.back()->setTimestamp(capture_time.toSec());
  {
    pcl::ScopeTime t("add frame occ time");
    map_ = map_ndt_->createOccupancyGrid(resolution_ / RESOLUTION_DIVISOR);
    repareOccupancyGrid();
  }
  means_ = map_ndt_->getMeansTransformed();
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
  //  std::sort(grids_.begin(), grids_.end(),
  //            [](const NDTGrid2DPtr &a, const NDTGrid2DPtr &b) {
  //              return a->operator<(*b);
  //            });
  map_ndt_->clear();
  for (const NDTGrid2DPtr &grid : grids_) {
    addToMap(*grid);
  }
  {
    pcl::ScopeTime t("add frame occ time");
    map_ = map_ndt_->createOccupancyGrid(resolution_ / RESOLUTION_DIVISOR);
    repareOccupancyGrid();
  }
  means_ = map_ndt_->getMeansTransformed();
  map_recalc_time_ = calc_time;
}

template <typename CellType, typename PointType>
void NDTMapper<CellType, PointType>::addToMap(const Grid &frame)
{
  map_ndt_->mergeInTraced(frame, true, true);
}

template <typename CellType, typename PointType>
void NDTMapper<CellType, PointType>::repareOccupancyGrid()
{
  cv::Mat image = frameToMat(map_);
  cv::Mat dst;
  int kernel_size = 1;
  cv::Mat kernel = cv::getStructuringElement(
      cv::MORPH_ELLIPSE, cv::Size(2 * kernel_size + 1, 2 * kernel_size + 1),
      cv::Point(kernel_size, kernel_size));
  // Apply erosion or dilation on the image
  cv::morphologyEx(image, dst, cv::MorphTypes::MORPH_CLOSE, kernel,
                   cv::Point(-1, -1), 2);
  // cv::erode(dst, image, kernel);
  map_.cells_ = matToFrame(dst, map_).cells_;
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

template <typename CellType, typename PointType>
cv::Mat
NDTMapper<CellType, PointType>::frameToMat(const OccupancyGrid &occ_grid) const
{
  cv::Mat frame_mat(occ_grid.height_, occ_grid.width_, CV_8UC1);
  for (size_t i = 0; i < occ_grid.cells_.size(); ++i) {
    if (occ_grid.cells_[i] == -1 || occ_grid.cells_[i] > 100)
      frame_mat.data[i] = 0;
    else if (occ_grid.cells_[i] < 50)
      frame_mat.data[i] = 0;
    else
      frame_mat.data[i] = occ_grid.cells_[i] + 150;
  }
  return frame_mat;
}

template <typename CellType, typename PointType>
OccupancyGrid NDTMapper<CellType, PointType>::matToFrame(
    const cv::Mat &img, const OccupancyGrid &original_grid) const
{
  OccupancyGrid grid;
  size_t idx = 0;
  for (int i = 0; i < img.rows; i++) {
    for (int j = 0; j < img.cols; j++) {
      auto pt = img.at<unsigned char>(i, j);
      if (original_grid.cells_[idx] == -1)
        grid.cells_.push_back(-1);
      else if (pt == 0)
        grid.cells_.push_back(0);
      else
        grid.cells_.push_back(pt - 150);
      ++idx;
    }
  }
  return grid;
}

// OccupancyGrid occ_grid = frame.createOccupancyGrid();
//   resolution_ = occ_grid.resolution_;
//    cv::Mat src = frameToMat(occ_grid);
//    // cv::Mat trans_src(frame_mat.size(),frame_mat.type());
//   float angle = eigt::getAngleDiffrence(Pose(0, 0, 0), frame->getOrigin());
//   float angle = eigt::getAngleDiffrence(Pose(0, 0, 0), frame.getOrigin());
//    // get rotation matrix for rotating the image around its center
//    cv::Point2f center(src.cols / 2.0f, src.rows / 2.0f);
//    cv::Mat rot = cv::getRotationMatrix2D(center, angle, 1.0);
//   // determine bounding rectangle
//   cv::Rect bbox = cv::RotatedRect(center, src.size(), angle).boundingRect();
//   // adjust transformation matrix
//   rot.at<double>(0, 2) += bbox.width / 2.0 - center.x;
//   rot.at<double>(1, 2) += bbox.height / 2.0 - center.y;
//    cv::Mat rotated;
//    cv::warpAffine(src, rotated, rot, bbox.size());
//    // blend rotated frame with global map
//   cv::Rect roi(occ_grid.origin_(0), occ_grid.origin_(1), rotated.cols,
//   auto coords_pair = calcCoordinates(occ_grid.origin_(0),
//   occ_grid.origin_(1),
//                                      width_, height_, occ_grid.resolution_);
//   cv::Rect roi(coords_pair.first, coords_pair.second, rotated.cols,
//                 rotated.rows);
//    addWeighted(map_(roi), 0.3, rotated, 0.7, 0.0, map_(roi));
//  }

}  // end of namespace slamuk
#endif
