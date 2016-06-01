#ifndef NDT_SCANMATCHING2D_ML_CORRELATIVE
#define NDT_SCANMATCHING2D_ML_CORRELATIVE

#include <ros/ros.h>
#include <pcl/registration/registration.h>
#include <pcl/filters/voxel_grid_covariance.h>
#include <Eigen/Dense>
#include <dynamic_slam_utils/eigen_tools.h>
#include <pcl/common/time.h>

namespace pcl
{
namespace ml_corr
{
struct SearchVoxel
{
  Eigen::Vector3f transform_;
  double score_;

  bool operator<(const SearchVoxel &rhs) const
  {
    return score_ < rhs.score_;
  }
};
////////////
struct IndexPoint
{
  int x_idx_;
  int y_idx_;
};
///////////////////////////////////////////LOOK UP TABLE UTILS /////////////////////////
class SmoothingKernel
{
public:
  typedef unsigned short CellType;

public:
  SmoothingKernel() : resolution_(0.5f), occupancy_(100)
  {
    initKernel(0.25f);
  }
  SmoothingKernel(float resolution, float std_deviation, short max_occupancy)
    : resolution_(resolution), occupancy_(max_occupancy)
  {
    initKernel(std_deviation);
  }
  int size() const
  {
    return size_;
  }
  int halfSize() const
  {
    return half_size_;
  }
  const CellType &operator[](size_t idx) const
  {
    return kernel_[idx];
  }
  /** return correct byte based on coordinates relative to the center of kernel
  * kernel:
  *  2    14   2           [-1,-1 ]    [-1,1]
  * 14   100  14 -------->        [0,0]
  *  2    14   2 [row,col] [1,-1]     [1,1]
  */
  const CellType & operator()(int row, int col) const{
    return kernel_[static_cast<size_t>((row + half_size_) * size_ + col+half_size_)];
  }

protected:
  float resolution_;
  short occupancy_;
  float std_deviation_;
  int size_;
  int half_size_;

  std::vector<CellType> kernel_;

  void initKernel(float std_deviation)
  {
    // check standard deviation limits
    std_deviation_ = std_deviation;
    if (std_deviation < 0.5f * resolution_)
      std_deviation_ = 0.5f * resolution_;
    if (std_deviation > 10 * resolution_)
      std_deviation_ = 10 * resolution_;
    // initialize size of kernel
    half_size_ =
        static_cast<int>(std::lround(2.0 * std_deviation_ / resolution_));
    size_ = 2 * half_size_ + 1;
    // fills kernel with normal distribution parameter
    kernel_.clear();
    for (int row = -half_size_; row <= half_size_; ++row) {
      for (int col = -half_size_; col <= half_size_; ++col) {
        double euclidean_dist =
            std::hypot(col * resolution_, row * resolution_);
        double val =
            std::exp(-0.5 * std::pow(euclidean_dist / std_deviation_, 2)) * occupancy_;
        assert(std::lround(val) < 255);
        kernel_.push_back(static_cast<CellType>(std::lround(val)));
      }
    }
  }
  friend std::ostream &operator<<(std::ostream &out, const SmoothingKernel &k);
};

std::ostream &operator<<(std::ostream &out, const SmoothingKernel &k)
{
  for (size_t i = 0; i < k.kernel_.size(); ++i) {
    if (i % k.size_ == 0)
      out << std::endl;
    out << std::setw(5)<< k.kernel_[i];
  }
  return out;
}

///////////////////////////////////////////
template <typename PointType>
class LookUpTable
{
  typedef short CellType;
public:
  LookUpTable()
    : cell_size_(0)
    , cell_size_half_(0)
    , minx_(0)
    , miny_(0)
    , maxx_(0)
    , maxy_(0)
    , cell_count_row_(0)
    , cell_count_col_(0)
    , border_size_ (0)
  {
  }
  void initGrid(const pcl::PointCloud<PointType> &target, float grid_step, float std_deviation);
  double getScore(const pcl::PointCloud<PointType> &pcl) const;
  double getScore(const std::vector<IndexPoint> &pcl) const;
  double getMaxScore()const;
  std::vector<IndexPoint> toIndexes(const pcl::PointCloud<PointType> &pcl) const;
  void moveIndexes(std::vector<IndexPoint> &indexes, int dx, int dy) const;
  void transformIndexes(const std::vector<IndexPoint> &source,
                        std::vector<IndexPoint> &out, float dx, float dy) const;

protected:
  float cell_size_;
  float cell_size_half_;
  // considered all points
  float minx_, miny_, maxx_, maxy_;
  size_t cell_count_row_;
  size_t cell_count_col_;
  size_t border_size_;
  size_t target_points_;
  const CellType OCCUPIED_CELL = 100;
  std::vector<CellType> table_;
  SmoothingKernel kernel_;

  void initGridWindowFce(const pcl::PointCloud<PointType> &target);
  size_t getCellIdx(const Eigen::Vector2d &pt) const;
  size_t getCellIdx(size_t row, size_t col) const;
  void applyKernel(size_t row, size_t col);
  // float calcCellScore(const std::vector<float> &grid, size_t i, size_t j) const;
  CellType getPtScore(const Eigen::Vector2d &pt) const;
  CellType getPtScore(const IndexPoint &pt) const;
  void getMinMax2D(const pcl::PointCloud<PointType> &pcl, float &minx,
                   float &miny, float &maxx, float &maxy) const;
  template <typename P>
  friend std::ostream &operator<<(std::ostream &out, const LookUpTable<P> &o);
};
////////////////////////////IMPLEMENTATION ///////////////////////
template <typename PointType>
void LookUpTable<PointType>::initGrid(const pcl::PointCloud<PointType> &target,
                                      float grid_step, float std_deviation)
{
  target_points_ = target.size();
  kernel_ =  SmoothingKernel(1/grid_step,std_deviation,OCCUPIED_CELL);
  border_size_ = kernel_.halfSize();
  cell_size_ = grid_step;
  cell_size_half_ = cell_size_ / 2;
  getMinMax2D(target, minx_, miny_, maxx_, maxy_);
  cell_count_row_ = std::ceil((maxx_ - minx_) / grid_step)+1+2*border_size_;
  cell_count_col_ = std::ceil((maxy_ - miny_) / grid_step)+1+2*border_size_;
  // kernel is used for smearing points in look up table

  std::cout<<"Smoothing kernel:\n"<<kernel_<<std::endl;
  ROS_DEBUG_STREAM("creating grid cells x: "<<cell_count_row_<< " y: "<<cell_count_col_);
  ROS_DEBUG_STREAM("values: minx: "<<minx_<<" maxx: "<<maxx_<<" miny "<<miny_<<" maxy "<<maxy_);
  initGridWindowFce(target);
  ROS_DEBUG_STREAM("grid initialized");

}

template <typename PointType>
double
LookUpTable<PointType>::getScore(const pcl::PointCloud<PointType> &pcl) const
{
  double res = 0;
  CellType temp = 0;
  size_t valid_pts = 0;
  for (size_t i = 0; i < pcl.size(); ++i) {
    temp = getPtScore(Eigen::Vector2d(pcl.points[i].x, pcl.points[i].y));
    if(temp == OCCUPIED_CELL)
      ++valid_pts;
    res += temp;
  }
  if(pcl.size() == 0)
    return 0;
  return res / (OCCUPIED_CELL *target_points_);
}

template <typename PointType>
double LookUpTable<PointType>::getScore(const std::vector<IndexPoint> &pcl) const
{
  double res = 0;
  CellType temp = 0;
  size_t valid_pts = 0;
  for (size_t i = 0; i < pcl.size(); ++i) {
    if (!(pcl[i].x_idx_ < 0 || pcl[i].y_idx_ < 0 || pcl[i].x_idx_ >= cell_count_row_ ||
      pcl[i].y_idx_ >= cell_count_col_)) {
      ++valid_pts;
    }
    temp = getPtScore(pcl[i]);
    //if(temp == OCCUPIED_CELL)
    //  ++valid_pts;
    res += temp;
  }
  if(pcl.size() == 0)
    return 0;
  //ROS_DEBUG_STREAM("valid points: "<<valid_pts <<" score: "<<res <<" (pts*score)/100 "<<res*valid_pts / 100);
  return res / (target_points_ * OCCUPIED_CELL);
}
template <typename PointType>
double LookUpTable<PointType>::getMaxScore()const{
  double sum = 0;
  size_t pts = 0;
  for(auto && cell : table_){
    sum+=cell;
    if(cell > 0)
      ++pts;
  }
  return sum ;/// (pts * OCCUPIED_CELL);
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
void LookUpTable<PointType>::initGridWindowFce(const pcl::PointCloud<PointType> &target)
{
  table_.clear();
  table_.reserve(cell_count_row_ * cell_count_col_);
  for (size_t i = 0; i < cell_count_row_ * cell_count_col_; ++i) {
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
      // if(getCellIdx(pt) > grid.size() -1){
      //   std::cout<<getCellIdx(pt)<<std::endl;
      //   std::cout<<pt.transpose()<<std::endl;
      // }
      size_t idx = getCellIdx(pt);
      table_[idx] =  OCCUPIED_CELL;
    }
  }
  //std::cout <<*this <<std::endl;
  ROS_DEBUG_STREAM("points used");
  // apply smoothing kernel
  for (size_t row = 1; row < cell_count_col_ - 1; ++row) {
    for (size_t col = 1; col < cell_count_row_ - 1; ++col) {
      size_t idx = getCellIdx(row,col);
      if(table_[idx] == OCCUPIED_CELL)
        applyKernel(row,col);

    }
  }
}

template <typename PointType>
size_t LookUpTable<PointType>::getCellIdx(const Eigen::Vector2d &pt) const
{
  size_t x_id = static_cast<size_t>(
      std::floor((pt(0) - minx_ + cell_size_half_) / cell_size_) + border_size_);
  size_t y_id = static_cast<size_t>(
      std::floor((-pt(1) + maxy_ + cell_size_half_) / cell_size_) + border_size_);
  return y_id * cell_count_row_ + x_id;
}

template <typename PointType>
size_t LookUpTable<PointType>::getCellIdx(size_t row, size_t col) const
{
  return row * cell_count_row_ + col;
}

template<typename PointType>
void LookUpTable<PointType>::applyKernel(size_t row, size_t col){
   // avoid applying kernel to the field at the edge - it is not necessary
   //size_t kernel_half_size = static_cast<size_t>(kernel_.halfSize());
   int half_kernel = kernel_.halfSize();
   // if (row < kernel_half_size || col < kernel_half_size ||
   //     col >= cell_count_row_ - 1 - kernel_half_size ||
   //     row >= cell_count_col_ - 1 - kernel_half_size)
   //   return;
   for(int r = - half_kernel; r <= half_kernel;++r){
    for(int c = -half_kernel; c <= half_kernel;++c){
      //
      size_t idx = getCellIdx(row+r,col+c);
      CellType kernel_val = kernel_(r,c);
      if(table_[idx] < kernel_val){
        table_[idx] =  kernel_val;
      }
    }
  }

}
// template <typename PointType>
// float LookUpTable<PointType>::calcCellScore(const std::vector<float> &grid,
//                                             size_t i, size_t j) const
// {
//   if (j == 0 || i == 0 || j >= cell_count_row_ - 2 || i >= cell_count_col_ - 2)
//     return 0;
//   if(grid[getCellIdx(i, j)] == 0)
//     return 0;
//   // getCellIdx(y,x)
//   // using gaussian kernel 3x3 for smoothing
//   float score = grid[getCellIdx(i - 1, j - 1)] * 0.075 +
//          grid[getCellIdx(i - 1, j)] * 0.124 +
//          grid[getCellIdx(i - 1, j + 1)] * 0.075 +
//          grid[getCellIdx(i, j - 1)] * 0.124 + grid[getCellIdx(i, j)] * 0.204 +
//          grid[getCellIdx(i, j + 1)] * 0.124 +
//          grid[getCellIdx(i + 1, j - 1)] * 0.075 +
//          grid[getCellIdx(i + 1, j)] * 0.124 +
//          grid[getCellIdx(i + 1, j + 1)] * 0.075;
//   //score = grid[getCellIdx(i, j)];
//   return score;
// }

template <typename PointType>
typename LookUpTable<PointType>::CellType
LookUpTable<PointType>::getPtScore(const Eigen::Vector2d &pt) const
{
  if (pt(0) < maxx_ && pt(1) < maxy_ && pt(0) > minx_ && pt(1) > miny_) {
    return table_[getCellIdx(pt)];
  }
  return 0;
}

template <typename PointType>
typename LookUpTable<PointType>::CellType
LookUpTable<PointType>::getPtScore(const IndexPoint &pt) const
{
  if (pt.x_idx_ < 0 || pt.y_idx_ < 0 || pt.x_idx_ >= cell_count_row_ ||
      pt.y_idx_ >= cell_count_col_) {
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
    if (i % o.cell_count_row_ == 0)
      out << std::endl;
    out << std::setw(5)<< o.table_[i];
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
  Eigen::MatrixXf rotation(4,4);
  rotation << co, -si,0,0,si,co,0,0,0,0,0,0,0,0,0,0;
  Eigen::MatrixXf points(4,source.size());
  outp.clear();
  outp.reserve(source.size());
  for (size_t i = 0; i < source.size(); ++i) {
    points.block<4,1>(0,i) << source[i].x, source[i].y,0,0;
     //  PointType pt = source[i];
     // pt.x = source[i].x * co - source[i].y * si;
     // pt.y = source[i].x * si + source[i].y * co;
     // outp.push_back(pt);
  }

  points =rotation * points;
  for (size_t i = 0; i < source.size(); ++i) {
    outp.push_back(PointType(points(0,i),points(1,i),0));
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
  : coarse_step_(0.25f)
  , fine_step_(0.1f)
  , coarse_rot_step_(0.2f)
  , fine_rot_step_(0.1f)
  , translation_range_(4.5f)
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
  transformPointCloud(*input_, guess_pcl, guess);
  // initialize lookup tables
  coarse_lookup_.initGrid(*target_, coarse_step_,0.5);
  fine_lookup_.initGrid(*target_, fine_step_,0.5);
  ROS_DEBUG_STREAM("grids initialized");
  //////////////////
   ml_corr::SearchVoxel voxel2;
   voxel2.transform_ << 0, 0, 0;
   voxel2.score_ = coarse_lookup_.getScore(guess_pcl);
   ROS_DEBUG_STREAM("[MultiLayerCorrelativeMatcher]: Coarse Voxel guess = trans: "
                    << voxel2.transform_.transpose()
                    << " score: " << voxel2.score_);

  voxel2.transform_ << 0, 0, 0;
  voxel2.score_ = coarse_lookup_.getScore(*target_);
  ROS_DEBUG_STREAM("[MultiLayerCorrelativeMatcher]: Coarse Voxel target= trans: "
                   << voxel2.transform_.transpose()
                   << " score: " << voxel2.score_);
  ROS_DEBUG_STREAM("[MultiLayerCorrelativeMatcher]: maximum on target: "<<coarse_lookup_.getMaxScore());
  // iterate over coarse grid /////////////////
  std::vector<ml_corr::SearchVoxel> search_voxels;
  size_t elements = static_cast<size_t>(
      std::ceil(std::pow((2 * translation_range_) / coarse_step_, 2) *
                (rotation_range_ / coarse_rot_step_)));
  size_t range_elements =
      static_cast<size_t>(std::floor(translation_range_ / coarse_step_));
  search_voxels.reserve(elements);

  std::vector<ml_corr::SearchVoxel> search_voxels_thread[4];
  // prepare all rotations
  std::vector<float> rotations;
  for(float theta = -rotation_range_; theta < rotation_range_;
         theta += coarse_rot_step_){
    rotations.push_back(theta);
  }
  {
    pcl::ScopeTime t_init("try total:");
//#pragma omp parallel num_threads(4)
 // {
   // pcl::ScopeTime t_init("try all translates:");
   pcl::PointCloud<PointSource> rot_pcl;
   std::vector<ml_corr::IndexPoint> index_pcl;
   std::vector<ml_corr::IndexPoint> transl_index_pcl;
#pragma omp parallel for private(rot_pcl,index_pcl, transl_index_pcl) num_threads(4)
    for (size_t theta_id = 0; theta_id < rotations.size(); ++theta_id) {
      int thread_id = omp_get_thread_num();
      // rotate point cloud
      ml_corr::rotatePointCloud(guess_pcl, rot_pcl, rotations[theta_id]);
      // float theta = 0;
      index_pcl = coarse_lookup_.toIndexes(rot_pcl);
      pcl::ScopeTime t_init("try single translates:");
      // iterate over all possible translation x values
      for (float dx = -translation_range_; dx < translation_range_;
           dx += coarse_step_) {
        // iterate over all possible y translation values
        for (float dy = -translation_range_; dy < translation_range_;
             dy += coarse_step_) {
          coarse_lookup_.transformIndexes(index_pcl, transl_index_pcl, dx, dy);
          ml_corr::SearchVoxel voxel;
          voxel.transform_ << dx, dy, rotations[theta_id];
          voxel.score_ = coarse_lookup_.getScore(transl_index_pcl);
          search_voxels_thread[thread_id].push_back(voxel);
        }
      }
    }
  //}  // end of parallel world
  std::cout<<"calc done"<<std::endl;
  // merge parallel results
  for (size_t i = 0; i < 4; ++i) {
    std::copy(search_voxels_thread[i].begin(), search_voxels_thread[i].end(),
              std::back_inserter(search_voxels));
  }
 }
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
  double best_score = 0;
  Eigen::Vector3f best_trans = search_voxels.back().transform_;
  Eigen::Matrix3f K = Eigen::Matrix3f::Identity();
  Eigen::Vector3f u = Eigen::Vector3f::Identity();
  float s = 0;
  ROS_DEBUG_STREAM("number of search voxels: "<<search_voxels.size());
  for (size_t i = search_voxels.size() - 1; i > search_voxels.size() - 20;
       --i) {
    ROS_DEBUG_STREAM("[MultiLayerCorrelativeMatcher]: Coarse Voxel = trans: "
                     << search_voxels[i].transform_.transpose()
                     << " score: " << search_voxels[i].score_);
  }

  //for(size_t i = 0; i < 8;++i){
  // while(true){
  //   best_voxel = search_voxels.back();
  //   search_voxels.pop_back();
  //   if(best_score > best_voxel.score_)
  //     break;

  //   for (float theta = best_voxel.transform_(2); theta < best_voxel.transform_(2) +coarse_rot_step_;
  //        theta += fine_rot_step_) {
  //     // rotate point cloud
  //     {
  //     pcl::ScopeTime t_init ("rotate pcl");
  //     ml_corr::rotatePointCloud(guess_pcl, rot_pcl, theta);
  //     //float theta = 0;
  //     index_pcl = fine_lookup_.toIndexes(rot_pcl);
  //     }
  //     {
  //       pcl::ScopeTime t_init ("try all translates:");
  //     // iterate over all possible translation x values
  //     for (float dx = best_voxel.transform_(0); dx < best_voxel.transform_(0)+fine_step_;
  //          dx += fine_step_) {
  //       // iterate over all possible y translation values
  //       for (float dy = best_voxel.transform_(1); dy < best_voxel.transform_(1)+fine_step_;
  //            dy += fine_step_) {
  //           fine_lookup_.transformIndexes(index_pcl, transl_index_pcl, dx, dy);
  //           double score = fine_lookup_.getScore(transl_index_pcl);
  //           if(score > best_score){
  //             best_score = score;
  //             best_trans<<dx,dy,theta;
  //           }
  //         }
  //       }
  //     }
  //   }
  //   ROS_DEBUG_STREAM("next fine iteration "<<best_score<<" ~ "<<best_voxel.score_);
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
