#ifndef NDT_GRID2D_NDT_MAPPER
#define NDT_GRID2D_NDT_MAPPER

#include <map>
#include <memory>
#include <Eigen/Dense>
#include <opencv/cv.h>
#include <dynamic_slam_utils/eigen_tools.h>
//#include <ndt_grid2d/types.h>

/*
FrameType needs to have implemented:
operator<
setOrigin(Pose)
setTimestamp(double)
getTimestamp(double)

*/
namespace slamuk
{

template <typename FrameType>
class NDTMapper
{
  typedef std::unique_ptr<FrameType> FrameTypePtr;
  typedef typename FrameType::Pose Pose;
  struct Comparer
  {
    bool operator()(const FrameTypePtr &a, const FrameTypePtr &b){
      return a->operator<(*b);
    }
  };

public:
  FrameType* addFrame(const FrameType &frame);
  FrameType* addFrame(FrameType && frame);
  void removeFrame(const FrameType & frame);
  void updateFrame(const FrameType & frame, const Pose &new_pose);
  FrameType* combineFrames(const FrameType & a, const FrameType & b);
  void recalc();

  const cv::Mat &getOccupancyMap() const;

protected:
  cv::Mat map_;
  size_t erased_elements = 0;
  std::vector<std::unique_ptr<FrameType>> grids_;

  // return id of frame if found. Otherwise exception
  size_t getFrameId(const FrameType & frame);
  void addToMap(FrameType * frame);
};

// ////////IMPLEMENTATION////
template <typename FrameType>
FrameType* NDTMapper<FrameType>::addFrame(const FrameType &frame)
{
  grids_.push_back(frame);
  addToMap(grids_.back());
  return grids_.back().get();
}

template <typename FrameType>
FrameType* NDTMapper<FrameType>::addFrame(FrameType &&frame)
{
  grids_.push_back(std::move(frame));
  addToMap(grids_.back());
  return grids_.back().get();
}

template <typename FrameType>
void NDTMapper<FrameType>::removeFrame(const FrameType &frame)
{

  grids_.erase();
}

template <typename FrameType>
void NDTMapper<FrameType>::updateFrame(const FrameType &frame, const Pose &pose)
{
  getFramePtr(id)->setOrigin(pose);
}

template <typename FrameType>
size_t NDTMapper<FrameType>::combineFrames(size_t first_id, size_t second_id){
  FrameType * frame_a = getFramePtr(first_id);
  FrameType * frame_b = getFramePtr(second_id);
  frame_a->mergeInTraced(*frame_b,true,true);
  frame_a->setTimestamp(frame_b -> getTimestamp());
  removeFrame(second_id);
}

template <typename FrameType>
void NDTMapper<FrameType>::recalc(){

}

template <typename FrameType>
void NDTMapper<FrameType>::addToMap(const FrameType & frame){
  OccupancyGrid occ_grid = frame.getOccupancyMap();
  cv::Mat src = frameToMat(frame);
  //cv::Mat trans_src(frame_mat.size(),frame_mat.type());
  float angle = eigt::getAngleDiffrence(Pose(0,0,0),frame->getOrigin());
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
  addWeighted( map_(roi), 0.3, rotated, 0.7, 0.0, map_(roi));
}

//// PROTECTED ////////
template <typename FrameType>
cv::Mat NDTMapper<FrameType>::frameToMat(const OccupancyGrid & occ_grid)
{
  cv::Mat frame_mat(occ_grid.height_, occ_grid.width_,CV_8UC1);
  for(size_t i = 0;i<occ_grid.cells_.size();++i){
    if(occ_grid.cells_[i] == -1 || occ_grid.cells_[i] > 100)
      frame_mat.data[i] = 205;
    else
      frame_mat.data[i] = occ_grid.cells_[i];
  }
}

template <typename FrameType>
size_t NDTMapper<FrameType>::getFrameId(const FrameType & frame)
{
  auto iter = grids_.find(id);
  if(iter == grids_.end())
    throw std::invalid_argument("Id of frame doesn't exist");
  return iter->get();
}


}// end of namespace slamuk
#endif
