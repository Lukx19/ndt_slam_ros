#ifndef NDT_SCANMATCHING2D_ML_CORRELATIVE
#define NDT_SCANMATCHING2D_ML_CORRELATIVE

#include <ros/ros.h>
#include <pcl/registration/registration.h>
#include <pcl/filters/voxel_grid_covariance.h>
#include <Eigen/Dense>

namespace pcl
{
namespace ml_corr
{
struct SearchVoxel
{
  Eigen::Vector3f transform_;
  float score_;

  bool operator<(const SearchVoxel &rhs) const
  {
    return score_ < rhs.score_;
  }
};

struct IndexPoint
{
  int x_idx_;
  int y_idx_;
};

template <typename PointType>
class LookUpTable
{
public:
  LookUpTable()
    : cell_size_(0)
    , cell_size_half_(0)
    , minx_(0)
    , miny_(0)
    , maxx_(0)
    , maxy_(0)
    , cell_count_x_(0)
    , cell_count_y_(0)
  {
  }
  void initGrid(const pcl::PointCloud<PointType> &target, float grid_step);
  float getScore(const pcl::PointCloud<PointType> &pcl) const;
  float getScore(const std::vector<IndexPoint> &pcl) const;
  std::vector<IndexPoint> toIndexes(const pcl::PointCloud<PointType> &pcl) const;
  void moveIndexes(std::vector<IndexPoint> &indexes, int dx, int dy) const;
  void transformIndexes(const std::vector<IndexPoint> &source,
                        std::vector<IndexPoint> &out, float dx, float dy) const;

protected:
  float cell_size_;
  float cell_size_half_;
  // considered all points
  float minx_, miny_, maxx_, maxy_;
  size_t cell_count_x_;
  size_t cell_count_y_;

  std::vector<float> table_;
  void initGridDistFce(const pcl::PointCloud<PointType> &target);
  void initGridWindowFce(const pcl::PointCloud<PointType> &target);
  size_t getCellIdx(const Eigen::Vector2d &pt) const;
  size_t getCellIdx(size_t idx_y, size_t idx_x) const;
  float calcCellScore(const std::vector<float> &grid, size_t i, size_t j) const;
  float getPtScore(const Eigen::Vector2d &pt) const;
  float getPtScore(const IndexPoint &pt) const;
  void getMinMax2D(const pcl::PointCloud<PointType> &pcl, float &minx,
                   float &miny, float &maxx, float &maxy) const;
  template <typename P>
  friend std::ostream &operator<<(std::ostream &out, const LookUpTable<P> &o);
};

template <typename PointType>
void LookUpTable<PointType>::initGrid(const pcl::PointCloud<PointType> &target,
                                      float grid_step)
{
  cell_size_ = grid_step;
  cell_size_half_ = cell_size_ / 2;
  getMinMax2D(target, minx_, miny_, maxx_, maxy_);
  cell_count_x_ = std::ceil((maxx_ - minx_) / grid_step);
  cell_count_y_ = std::ceil((maxy_ - miny_) / grid_step);
  ROS_DEBUG_STREAM("creating grid"<<cell_count_x_<< " y: "<<cell_count_y_);
  initGridWindowFce(target);
  //initGridDistFce(target);
  ROS_DEBUG_STREAM("grid initialized");

}

template <typename PointType>
float
LookUpTable<PointType>::getScore(const pcl::PointCloud<PointType> &pcl) const
{
  float res = 0;
  for (size_t i = 0; i < pcl.points.size(); ++i) {
    res += getPtScore(Eigen::Vector2d(pcl.points[i].x, pcl.points[i].y));
  }
  return res;
}

template <typename PointType>
float LookUpTable<PointType>::getScore(const std::vector<IndexPoint> &pcl) const
{
  float res = 0;
  float temp = 0;
  // size_t valid_pts = 0;
  for (size_t i = 0; i < pcl.size(); ++i) {
    temp = getPtScore(pcl[i]);
    // if(temp > 0)
    //  ++valid_pts;
    res += temp;
  }
  // ROS_DEBUG_STREAM("valid points: "<<valid_pts);
  return res;
}
template <typename PointType>
std::vector<IndexPoint>
LookUpTable<PointType>::toIndexes(const pcl::PointCloud<PointType> &pcl) const
{
  std::vector<IndexPoint> indexes;
  for (size_t i = 0; i < pcl.points.size(); ++i) {
    IndexPoint pt;
    pt.x_idx_ = static_cast<int>(
        std::floor((pcl.points[i].x - minx_ + cell_size_half_) / cell_size_));
    pt.y_idx_ = static_cast<int>(
        std::floor((-pcl.points[i].y + maxy_ + cell_size_half_) / cell_size_));
    indexes.push_back(pt);
  }
  return indexes;
}

template <typename PointType>
void LookUpTable<PointType>::moveIndexes(std::vector<IndexPoint> &indexes,
                                         int dx, int dy) const
{
  for (size_t i = 0; i < indexes.size(); ++i) {
    indexes[i].x_idx_ += dx;
    indexes[i].y_idx_ += dy;
  }
}

template <typename PointType>
void
LookUpTable<PointType>::transformIndexes(const std::vector<IndexPoint> &source,
                                         std::vector<IndexPoint> &out, float dx,
                                         float dy) const
{
  out.clear();
  out.reserve(source.size());
  int add_x = static_cast<int>(std::floor(dx / cell_size_));
  int add_y = -static_cast<int>(std::floor(dy / cell_size_));
  // std::cout<<"transformPCLIndexes "<<add_x<<" y: "<<add_y<<"\n";
  for (size_t i = 0; i < source.size(); ++i) {
    IndexPoint new_pt;
    new_pt.x_idx_ = source[i].x_idx_ + add_x;
    new_pt.y_idx_ = source[i].y_idx_ + add_y;
    out.push_back(new_pt);
  }
}


/////////////////////////////////// PROTECTED ///////////////////////////////
template <typename PointType>
void LookUpTable<PointType>::initGridDistFce(const pcl::PointCloud<PointType> &target)
{
  std::vector<float> grid;
  grid.reserve(cell_count_x_ * cell_count_y_);
  table_.clear();
  table_.reserve(cell_count_x_ * cell_count_y_);
  for (size_t i = 0; i < cell_count_x_ * cell_count_y_; ++i) {
    grid.push_back(0);
    table_.push_back(0);
  }
  ROS_DEBUG_STREAM("structures created");
  // initialize target grid
  // every point is projected to certain cell. If is cell in the grid we set
  // this cell 1
  float sigma  = 0.01f;
  float norm_c = 1 / (std::sqrt(2*3.14159265359f) * sigma);
  for (size_t i = 0; i < target.size(); ++i) {
    if (target[i].x < maxx_ && target[i].y < maxy_ &&
        target[i].x > minx_ && target[i].y > miny_) {
      size_t idx = getCellIdx(Eigen::Vector2d(target[i].x, target[i].y));
      float dist = std::sqrt(std::pow(target[i].x,2) + std::pow(target[i].y,2));
      float weight = std::exp(-0.5f * std::pow((dist/sigma),2.0f));// * norm_c;
       ROS_DEBUG_STREAM("weight"<<weight<<" norm "<<norm_c<< " dist "<<dist);
      grid[idx] +=0.01;
    }
  }
  ROS_DEBUG_STREAM("points used");
  for (size_t i = 0; i < target.size(); ++i) {
    size_t idx = getCellIdx(Eigen::Vector2d(target[i].x, target[i].y));
    for(size_t grid_id = 0; grid_id < grid.size();++grid_id){
      //table_[grid_id] = std::min(grid[idx] + table_[grid_id],39.894f);
      table_[idx]+=grid[grid_id];
    }
  }
}

template <typename PointType>
void LookUpTable<PointType>::initGridWindowFce(const pcl::PointCloud<PointType> &target)
{
  std::vector<float> grid;
  grid.reserve(cell_count_x_ * cell_count_y_);
  table_.clear();
  table_.reserve(cell_count_x_ * cell_count_y_);
  for (size_t i = 0; i < cell_count_x_ * cell_count_y_; ++i) {
    grid.push_back(0);
    table_.push_back(0);
  }
  ROS_DEBUG_STREAM("structures created");
  // initialize target grid
  // every point is projected to certain cell. If is cell in the grid we set
  // this cell 1
  for (size_t i = 0; i < target.size(); ++i) {
    if (target[i].x < maxx_ && target[i].y < maxy_ &&
        target[i].x > minx_ && target[i].y > miny_) {
      Eigen::Vector2d pt(target[i].x, target[i].y);
      std::cout<<getCellIdx(pt)<<std::endl;
      //size_t idx = getCellIdx(pt);
      grid[0] =  1;
    }
  }
  ROS_DEBUG_STREAM("points used");
  // initialize gaussian values
  // for (size_t i = 1; i < cell_count_y_ - 1; ++i) {
  //   for (size_t j = 1; j < cell_count_x_ - 1; ++j) {
  //     //ROS_DEBUG_STREAM("pt_pred:"<<i<<"~"<<j);
  //     table_[i * cell_count_x_ + j] = calcCellScore(grid, i, j);
  //     //ROS_DEBUG_STREAM("pt_po:"<<i<<"~"<<j);
  //   }
  // }
}

template <typename PointType>
size_t LookUpTable<PointType>::getCellIdx(const Eigen::Vector2d &pt) const
{
  size_t x_id = static_cast<size_t>(
      std::floor((pt(0) - minx_ + cell_size_half_) / cell_size_));
  size_t y_id = static_cast<size_t>(
      std::floor((-pt(1) + maxy_ + cell_size_half_) / cell_size_));
  return y_id * cell_count_x_ + x_id;
}

template <typename PointType>
size_t LookUpTable<PointType>::getCellIdx(size_t idx_y, size_t idx_x) const
{
  return idx_y * cell_count_x_ + idx_x;
}

template <typename PointType>
float LookUpTable<PointType>::calcCellScore(const std::vector<float> &grid,
                                            size_t i, size_t j) const
{
  if (j == 0 || i == 0 || j >= cell_count_x_ - 2 || i >= cell_count_y_ - 2)
    return 0;
  // if(grid[getCellIdx(i, j)] == 0)
  //  return 0;
  // getCellIdx(y,x)
  float score = grid[getCellIdx(i - 1, j - 1)] * 0.075 +
         grid[getCellIdx(i - 1, j)] * 0.124 +
         grid[getCellIdx(i - 1, j + 1)] * 0.075 +
         grid[getCellIdx(i, j - 1)] * 0.124 + grid[getCellIdx(i, j)] * 0.204 +
         grid[getCellIdx(i, j + 1)] * 0.124 +
         grid[getCellIdx(i + 1, j - 1)] * 0.075 +
         grid[getCellIdx(i + 1, j)] * 0.124 +
         grid[getCellIdx(i + 1, j + 1)] * 0.075;
  //score = grid[getCellIdx(i, j)];
  return score;
}

template <typename PointType>
float LookUpTable<PointType>::getPtScore(const Eigen::Vector2d &pt) const
{
  if (pt(0) < maxx_ && pt(1) < maxy_ && pt(0) > minx_ && pt(1) > miny_) {
    return table_[getCellIdx(pt)];
  }
  return 0;
}

template <typename PointType>
float LookUpTable<PointType>::getPtScore(const IndexPoint &pt) const
{
  if (pt.x_idx_ < 0 || pt.y_idx_ < 0 || pt.x_idx_ >= cell_count_x_ ||
      pt.y_idx_ >= cell_count_y_) {
    return 0;
  }
  return table_[getCellIdx(pt.y_idx_, pt.x_idx_)];
}

template <typename PointType>
void LookUpTable<PointType>::getMinMax2D(const pcl::PointCloud<PointType> &pcl,
                                         float &minx, float &miny, float &maxx,
                                         float &maxy) const
{
  minx = miny = maxx = maxy = 0;
  for (size_t i = 0; i < pcl.points.size(); ++i) {
    if (pcl.points[i].x < minx)
      minx = pcl.points[i].x;
    if (pcl.points[i].x > maxx)
      maxx = pcl.points[i].x;
    if (pcl.points[i].y < miny)
      miny = pcl.points[i].y;
    if (pcl.points[i].y > maxy)
      maxy = pcl.points[i].y;
  }
}

template <typename P>
std::ostream &operator<<(std::ostream &out, const LookUpTable<P> &o)
{
  for (size_t i = 0; i < o.table_.size(); ++i) {
    if (i % o.cell_count_x_ == 0)
      out << std::endl;
    out << " " << static_cast<float>(o.table_[i]);
  }
  return out;
}

template <typename PointType>
void rotatePointCloud(const pcl::PointCloud<PointType> &source,
                      pcl::PointCloud<PointType> &outp, float angle)
{
  // x' = x cos f - y sin f
  // y' = y cos f + x sin f
  float si = std::sin(angle);
  float co = std::cos(angle);
  outp.clear();
  outp.reserve(source.size());
  for (size_t i = 0; i < source.size(); ++i) {
    PointType pt = source[i];
    pt.x = source[i].x * co - source[i].y * si;
    pt.y = source[i].y * co - source[i].x * si;
    outp.push_back(pt);
  }
}

template <typename PointType>
void translatePointCloud(const pcl::PointCloud<PointType> &source,
                         pcl::PointCloud<PointType> &outp, float dx, float dy)
{
  outp.clear();
  outp.reserve(source.size());
  for (size_t i = 0; i < source.size(); ++i) {
    PointType pt = source[i];
    pt.x += dx;
    pt.y += dy;
    outp.push_back(pt);
  }
}

}  // end of namespace ml_corr

template <typename PointSource, typename PointTarget>
class MultiLayerCorrelativeMatcher
    : public Registration<PointSource, PointTarget>
{
protected:
  typedef typename Registration<PointSource, PointTarget>::PointCloudSource
      PclSource;
  typedef typename PclSource::Ptr PclSourcePtr;
  typedef typename PclSource::ConstPtr PclSourceConstPtr;

  typedef typename Registration<PointSource, PointTarget>::PointCloudTarget
      PclTarget;
  typedef typename PclTarget::Ptr PclTargetPtr;
  typedef typename PclTarget::ConstPtr PclTargetConstPtr;

  typedef PointIndices::Ptr PointIndicesPtr;
  typedef PointIndices::ConstPtr PointIndicesConstPtr;

  /** \brief Typename of searchable voxel grid containing mean and covariance.
   */
  typedef VoxelGridCovariance<PointTarget> TargetGrid;
  typedef VoxelGridCovariance<PointSource> SourceGrid;

  /** \brief Typename of const pointer to searchable voxel grid. */
  typedef const TargetGrid *TargetGridConstPtr;
  /** \brief Typename of const pointer to searchable voxel grid leaf. */
  typedef typename TargetGrid::LeafConstPtr TargetGridLeafConstPtr;

public:
  typedef boost::shared_ptr<
      MultiLayerCorrelativeMatcher<PointSource, PointTarget>> Ptr;
  typedef boost::shared_ptr<
      const MultiLayerCorrelativeMatcher<PointSource, PointTarget>> ConstPtr;
  typedef Eigen::Vector3d VectorTrans;

  MultiLayerCorrelativeMatcher();

  /** \brief Empty destructor */
  virtual ~MultiLayerCorrelativeMatcher()
  {
  }

  inline void setCoarseStep(float step)
  {
    coarse_step_ = step;
  }

  inline float getCoarseStep() const
  {
    return coarse_step_;
  }

  inline Eigen::Matrix3d getCovariance() const
  {
    return covariance_;
  }

  inline Eigen::Matrix3d getInformMatrix() const
  {
    return inform_matrix_;
  }

protected:
  using Registration<PointSource, PointTarget>::input_;
  using Registration<PointSource, PointTarget>::target_;
  using Registration<PointSource, PointTarget>::final_transformation_;
  using Registration<PointSource, PointTarget>::converged_;

  Eigen::Matrix3d covariance_;
  Eigen::Matrix3d inform_matrix_;

  float coarse_step_;
  float fine_step_;
  float coarse_rot_step_;
  float fine_rot_step_;
  float translation_range_;  // search from [-range,range]
  float rotation_range_;     // maximum is [-Pi,Pi]
  float MIN_OVERLAP_SCORE = 0;

  ml_corr::LookUpTable<PointTarget> coarse_lookup_;
  ml_corr::LookUpTable<PointTarget> fine_lookup_;

  /** \brief Estimate the transformation and returns the transformed source
   * (input) as output.
    * \param[out] output the resultant input transfomed point cloud dataset
    */
  virtual void computeTransformation(PclSource &output)
  {
    computeTransformation(output, Eigen::Matrix4f::Identity());
  }

  /** \brief Estimate the transformation and returns the transformed source
   * (input) as output.
    * \param[out] output the resultant input transfomed point cloud dataset
    * \param[in] guess the initial gross estimation of the transformation
    */
  virtual void computeTransformation(PclSource &output,
                                     const Eigen::Matrix4f &guess);

  template <typename T = float, typename In>
  Eigen::Matrix<T, 4, 4> vecToMat(const Eigen::Matrix<In, 3, 1> &trans) const;
};

template <typename PointSource, typename PointTarget>
MultiLayerCorrelativeMatcher<PointSource,
                             PointTarget>::MultiLayerCorrelativeMatcher()
  : coarse_step_(0.5f)
  , fine_step_(0.05f)
  , coarse_rot_step_(0.1f)
  , fine_rot_step_(0.01f)
  , translation_range_(2.5f)
  , rotation_range_(2.5f)
{
  converged_ = false;
}

template <typename PointSource, typename PointTarget>
void
MultiLayerCorrelativeMatcher<PointSource, PointTarget>::computeTransformation(
    PclSource &output, const Eigen::Matrix4f &guess)
{
  pcl::PointCloud<PointSource> guess_pcl;
  pcl::PointCloud<PointSource> rot_pcl;
  std::vector<ml_corr::IndexPoint> index_pcl;
  std::vector<ml_corr::IndexPoint> transl_index_pcl;
  transformPointCloud(*input_, guess_pcl, guess);
  // initialize lookup tables
  coarse_lookup_.initGrid(*target_, coarse_step_);
  //fine_lookup_.initGrid(*target_, fine_step_);
  ROS_DEBUG_STREAM("grids initialized");
  //////////////////
  // ml_corr::SearchVoxel voxel2;
  // voxel2.transform_ << 0, 0, 0;
  // voxel2.score_ = coarse_lookup_.getScore(guess_pcl);
  // ROS_DEBUG_STREAM("[MultiLayerCorrelativeMatcher]: Coarse Voxel guess = trans: "
  //                  << voxel2.transform_.transpose()
  //                  << " score: " << voxel2.score_);

  // voxel2.transform_ << 0, 0, 0;
  // voxel2.score_ = coarse_lookup_.getScore(coarse_lookup_.toIndexes(*target_));
  // ROS_DEBUG_STREAM("[MultiLayerCorrelativeMatcher]: Coarse Voxel target= trans: "
  //                  << voxel2.transform_.transpose()
  //                  << " score: " << voxel2.score_);
  // iterate over coarse grid /////////////////
  std::vector<ml_corr::SearchVoxel> search_voxels;
  size_t elements = static_cast<size_t>(
      std::ceil(std::pow((2 * translation_range_) / coarse_step_, 2) *
                (rotation_range_ / coarse_rot_step_)));
  size_t range_elements =
      static_cast<size_t>(std::floor(translation_range_ / coarse_step_));
  std::cout << "range elems " << range_elements << "\n";
  search_voxels.reserve(elements);
  // for (float theta = -rotation_range_; theta < rotation_range_;
  //      theta += coarse_rot_step_) {
  //   // rotate point cloud
  //   ml_corr::rotatePointCloud(guess_pcl, rot_pcl, theta);
  //   index_pcl = coarse_lookup_.toIndexes(rot_pcl);
  //   // iterate over all possible translation x values
  //   for (float dx = -translation_range_; dx < translation_range_;
  //        dx += coarse_step_) {
  //     // iterate over all possible y translation values
  //     for (float dy = -translation_range_; dy < translation_range_;
  //          dy += coarse_step_) {
  //       // ml_corr::translatePointCloud(*target_,rot_pcl,dx,dy);
  //       coarse_lookup_.transformIndexes(index_pcl, transl_index_pcl, dx, dy);
  //       // std::cout<<"x: "<<index_pcl[0].x_idx_<< "~"<<
  //       // transl_index_pcl[0].x_idx_<<"\n";
  //       // std::cout<<"y "<<index_pcl[0].y_idx_<< "~"<<
  //       // transl_index_pcl[0].y_idx_<<"\n";
  //       ml_corr::SearchVoxel voxel;
  //       voxel.transform_ << dx, dy, theta;
  //       voxel.score_ = coarse_lookup_.getScore(transl_index_pcl);
  //       // ROS_DEBUG_STREAM("[MultiLayerCorrelativeMatcher]: Coarse Voxel =
  //       // trans:"
  //       //                  << voxel.transform_.transpose()
  //       //                  << " score: " << voxel.score_);
  //       search_voxels.push_back(voxel);
  //       //++y_adds;
  //     }
  //     // coarse_lookup_.moveIndexes(index_pcl,0.0f,-y_adds);
  //   }
  // }

  // search_voxels.clear();

  //  std::vector<ml_corr::IndexPoint> transl_index_pcl;
  // for(float theta = -rotation_range_;
  // theta<rotation_range_;theta+=coarse_rot_step_){
  //   //rotate point cloud
  //    ml_corr::rotatePointCloud(guess_pcl,rot_pcl,theta);
  //    index_pcl = coarse_lookup_.toIndexes(rot_pcl);
  //    coarse_lookup_.transformIndexes(index_pcl,transl_index_pcl,-translation_range_,-translation_range_);
  //    ml_corr::SearchVoxel voxel;
  //    voxel.transform_<<0,0,theta;
  //    voxel.score_ = coarse_lookup_.getScore(transl_index_pcl);
  //    search_voxels.push_back(voxel);
  // }
  // std::cout<<search_voxels.size()<<"\n";

  // find voxel with best score = best chance for good match
  std::sort(search_voxels.begin(), search_voxels.end());
  ml_corr::SearchVoxel best_voxel;
  // test if this score is sufficent-> small score means little or no overlap
  if (search_voxels.back().score_ < MIN_OVERLAP_SCORE) {
    converged_ = false;
    return;
  }
  ROS_DEBUG_STREAM("[MultiLayerCorrelativeMatcher]: delta_trans:"
                   << search_voxels.back().transform_.transpose());
  float best_score = 0;
  Eigen::Vector3f best_trans = search_voxels.back().transform_;
  Eigen::Matrix3f K = Eigen::Matrix3f::Identity();
  Eigen::Vector3f u = Eigen::Vector3f::Identity();
  float s = 0;
  ROS_DEBUG_STREAM("number of search voxels: "<<search_voxels.size());
  // for (size_t i = search_voxels.size() - 1; i > search_voxels.size() - 20;
  //      --i) {
  //   ROS_DEBUG_STREAM("[MultiLayerCorrelativeMatcher]: Coarse Voxel = trans: "
  //                    << search_voxels[i].transform_.transpose()
  //                    << " score: " << search_voxels[i].score_);
  // }

  // for(size_t i = 0; i < 8;++i){
  //   best_voxel = search_voxels[i];
  //   //search_voxels.pop_back();
  //   if(best_score > best_voxel.score_)
  //     break;
  //   for (float theta = best_voxel.transform_(2);
  //        theta < best_voxel.transform_(2) + coarse_rot_step_;
  //        theta += fine_rot_step_){
  //     // rotate point cloud
  //     Eigen::Matrix4f trans_mat;
  //    Eigen::Rotation2D<float> rotate_mat(theta);
  //    trans_mat.setIdentity();
  //    trans_mat.block(0,0,2,2) = rotate_mat.toRotationMatrix();
  //    trans_mat(0,3) = trans_mat(1,3) = -translation_range_;
  //    transformPointCloud(guess_pcl, rot_pcl, trans_mat);
  //    index_pcl = fine_lookup_.toIndexes(rot_pcl);
  //     // iterate over all possible translation x values in voxel
  //     for (float dx = best_voxel.transform_(0);
  //          dx < best_voxel.transform_(0) + coarse_step_; dx += fine_step_) {
  //       //iterate over all possible y translation values in voxel
  //       fine_lookup_.moveIndexes(index_pcl,1,0.0f);
  //       size_t y_adds = 0;
  //       for(float dy = best_voxel.transform_(1);
  //          dy < best_voxel.transform_(1) + coarse_step_; dy += fine_step_){
  //         fine_lookup_.moveIndexes(index_pcl,0.0f,1);
  //         float score = fine_lookup_.getScore(index_pcl);
  //         Eigen::Vector3f xi;
  //         xi << dx,dy,theta;
  //         if(score > best_score){
  //           best_score = score;
  //           best_trans = xi;
  //         }
  //         // adding to vars for covarianve calculation
  //         K+=(xi*xi.transpose()) * score;
  //         u+=xi*score;
  //         s+=score;
  //          ++y_adds;
  //       }
  //       coarse_lookup_.moveIndexes(index_pcl,0.0f,-y_adds);
  //     }
  //   }
  //   ROS_DEBUG_STREAM("next fine iteration "<<best_score<<" ~
  //   "<<search_voxels[i].score_);
  //   //break;
  // }
  ROS_DEBUG_STREAM("[MultiLayerCorrelativeMatcher]: final_delta_trans:"
                   << best_trans.transpose());
  ROS_DEBUG_STREAM("\n"<<coarse_lookup_);
  // calculate covariance
  covariance_ =
      ((1 / s) * K - (1 / std::pow(s, 2) * (u * u.transpose()))).cast<double>();
  inform_matrix_ = covariance_.inverse();

  // output data
  converged_ = true;
  final_transformation_ = vecToMat(best_trans);  // * guess;
  transformPointCloud(*input_, output, final_transformation_);
}

template <typename PointSource, typename PointTarget>
template <typename T, typename In>
Eigen::Matrix<T, 4, 4>
MultiLayerCorrelativeMatcher<PointSource, PointTarget>::vecToMat(
    const Eigen::Matrix<In, 3, 1> &trans) const
{
  Eigen::Matrix<T, 4, 4> trans_mat = Eigen::Matrix<T, 4, 4>::Identity();

  trans_mat.block(0, 0, 3, 3).matrix() =
      Eigen::Matrix<T, 3, 3>(Eigen::AngleAxis<T>(
          static_cast<T>(trans(2)), Eigen::Matrix<T, 3, 1>::UnitZ()));

  trans_mat.block(0, 3, 3, 1).matrix() = Eigen::Matrix<T, 3, 1>(
      static_cast<T>(trans(0)), static_cast<T>(trans(1)), 0.0);

  return trans_mat;
}

}  // end of pcl namespace

#endif
