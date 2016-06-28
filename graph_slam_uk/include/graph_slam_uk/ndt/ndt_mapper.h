#ifndef GRAPH_SLAM_UK_NDT_MAPPER
#define GRAPH_SLAM_UK_NDT_MAPPER

#include <map>
#include <memory>
#include <Eigen/Dense>
#include <opencv/cv.h>
#include <graph_slam_uk/utils/eigen_tools.h>
#include <graph_slam_uk/ndt/ndt_grid2d_interface.h>

/*
INDTGrid2D needs to have implemented:
operator<
setOrigin(Pose)
setTimestamp(double)
getTimestamp(double)

*/
namespace slamuk
{
// template <typename INDTGrid2D>
class NDTMapper
{
  // typedef INDTGrid2D INDTGrid2D;
  typedef std::unique_ptr<INDTGrid2D> INDTGrid2DPtr;
  typedef typename INDTGrid2D::Pose Pose;
  struct Comparer
  {
    bool operator()(const INDTGrid2DPtr &a, const INDTGrid2DPtr &b)
    {
      return a->operator<(*b);
    }
  };

public:
  // INDTGrid2D *addFrame(const INDTGrid2D &frame);
  INDTGrid2D *addFrame(INDTGrid2DPtr &&frame);
  void removeFrame(const INDTGrid2D *frame);
  void updateFrame(const INDTGrid2D *frame, const Pose &new_pose);
  INDTGrid2D *combineFrames(const INDTGrid2D *a, const INDTGrid2D *b);
  void recalc();

  const cv::Mat &getOccupancyMap() const;

protected:
  cv::Mat map_;
  size_t erased_elements = 0;
  std::vector<INDTGrid2DPtr> grids_;

  // return id of frame if found. Otherwise exception
  size_t getFrameId(const INDTGrid2D &frame);
  void addToMap(INDTGrid2D *frame);
};

// ////////IMPLEMENTATION////
// template <typename INDTGrid2D>
// INDTGrid2D *NDTMapper<INDTGrid2D>::addFrame(const INDTGrid2D &frame)
// {
//   grids_.push_back(frame);
//   addToMap(grids_.back());
//   return grids_.back().get();
// }

INDTGrid2D *NDTMapper::addFrame(std::unique_ptr<INDTGrid2D> &&frame)
{
  grids_.push_back(std::move(frame));
  addToMap(grids_.back());
  return grids_.back().get();
}

void NDTMapper::removeFrame(const INDTGrid2D *frame)
{
}

void NDTMapper::updateFrame(const INDTGrid2D &frame, const Pose &pose)
{
  getFramePtr(id)->setOrigin(pose);
}

size_t NDTMapper::combineFrames(size_t first_id, size_t second_id)
{
  INDTGrid2D *frame_a = getFramePtr(first_id);
  INDTGrid2D *frame_b = getFramePtr(second_id);
  frame_a->mergeInTraced(*frame_b, true, true);
  frame_a->setTimestamp(frame_b->getTimestamp());
  removeFrame(second_id);
}

void NDTMapper::recalc()
{
}

void NDTMapper::addToMap(const INDTGrid2D *frame)
{
  OccupancyGrid occ_grid = frame.getOccupancyMap();
  cv::Mat src = frameToMat(occ_grid);
  // cv::Mat trans_src(frame_mat.size(),frame_mat.type());
  float angle = eigt::getAngleDiffrence(Pose(0, 0, 0), frame->getOrigin());
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
  cv::Rect roi(occ_grid.origin_(0), occ_grid.origin_(1), rotated.cols,
               rotated.rows);
  addWeighted(map_(roi), 0.3, rotated, 0.7, 0.0, map_(roi));
}

//// PROTECTED ////////

cv::Mat NDTMapper::frameToMat(const OccupancyGrid &occ_grid)
{
  cv::Mat frame_mat(occ_grid.height_, occ_grid.width_, CV_8UC1);
  for (size_t i = 0; i < occ_grid.cells_.size(); ++i) {
    if (occ_grid.cells_[i] == -1 || occ_grid.cells_[i] > 100)
      frame_mat.data[i] = 205;
    else
      frame_mat.data[i] = occ_grid.cells_[i];
  }
}

size_t NDTMapper::getFrameId(const INDTGrid2D *frame)
{
  auto iter = grids_.find(frame);
  if (iter == grids_.end())
    throw std::invalid_argument("Id of frame doesn't exist");
  return iter->get();
}

}  // end of namespace slamuk
#endif
