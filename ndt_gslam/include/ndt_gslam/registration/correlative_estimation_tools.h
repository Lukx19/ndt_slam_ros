#ifndef NDT_GSLAM_CORRELATIVE_ESTIMATION_TOOLS
#define NDT_GSLAM_CORRELATIVE_ESTIMATION_TOOLS

#include <math.h>
#include <pcl/point_cloud.h>
#include <ros/ros.h>
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

namespace pcl
{
namespace ml_corr
{
struct SearchVoxel {
  Eigen::Vector3f transform_;
  double score_;

  SearchVoxel() : transform_(0, 0, 0), score_(0)
  {
  }
  bool operator<(const SearchVoxel &rhs) const
  {
    return score_ < rhs.score_;
  }
};
////////////
struct IndexPoint {
  int x_idx_;
  int y_idx_;
};
///////////////////////////////////////////SMOOTHING KERNEL//////////////
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
  const CellType &operator()(int row, int col) const
  {
    return kernel_[static_cast<size_t>((row + half_size_) * size_ + col +
                                       half_size_)];
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
            std::exp(-0.5 * std::pow(euclidean_dist / std_deviation_, 2)) *
            occupancy_;
        assert(std::lround(val) < 255);
        kernel_.push_back(static_cast<CellType>(std::lround(val)));
      }
    }
  }
  friend std::ostream &operator<<(std::ostream &out, const SmoothingKernel &k);
};

inline std::ostream &operator<<(std::ostream &out, const SmoothingKernel &k)
{
  for (size_t i = 0; i < k.kernel_.size(); ++i) {
    if (i % k.size_ == 0)
      out << std::endl;
    out << std::setw(5) << k.kernel_[i];
  }
  return out;
}

////////////////////////////////////////////// LOOK UP TABLE///////
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
    , border_size_(0)
  {
  }
  void initGrid(const pcl::PointCloud<PointType> &target, float grid_step,
                float std_deviation);
  double getScore(const pcl::PointCloud<PointType> &pcl) const;
  double getScore(const std::vector<IndexPoint> &pcl) const;
  double getMaxScore() const;
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
  size_t getCellIdx(const Eigen::Vector2f &pt) const;
  size_t getCellIdx(size_t row, size_t col) const;
  void applyKernel(size_t row, size_t col);
  // float calcCellScore(const std::vector<float> &grid, size_t i, size_t j)
  // const;
  float getPtScore(const Eigen::Vector2f &pt) const;
  float getPtScore(const IndexPoint &pt) const;
  void getMinMax2D(const pcl::PointCloud<PointType> &pcl, float &minx,
                   float &miny, float &maxx, float &maxy) const;
  template <typename P>
  friend std::ostream &operator<<(std::ostream &out, const LookUpTable<P> &o);
};
////////////////////////////IMPLEMENTATION
template <typename PointType>
void LookUpTable<PointType>::initGrid(const pcl::PointCloud<PointType> &target,
                                      float grid_step, float std_deviation)
{
  target_points_ = target.size();
  kernel_ = SmoothingKernel(1 / grid_step, std_deviation, OCCUPIED_CELL);
  border_size_ = kernel_.halfSize();
  cell_size_ = grid_step;
  cell_size_half_ = cell_size_ / 2;
  getMinMax2D(target, minx_, miny_, maxx_, maxy_);
  cell_count_row_ =
      std::ceil((maxx_ - minx_) / grid_step) + 1 + 2 * border_size_;
  cell_count_col_ =
      std::ceil((maxy_ - miny_) / grid_step) + 1 + 2 * border_size_;
  // kernel is used for smearing points in look up table

  // std::cout << "Smoothing kernel:\n" << kernel_ << std::endl;
  ROS_DEBUG_STREAM("[CorrelativeEstimationGrid]:creating grid cells x: "
                   << cell_count_row_ << " y: " << cell_count_col_);
  initGridWindowFce(target);
}

template <typename PointType>
double
LookUpTable<PointType>::getScore(const pcl::PointCloud<PointType> &pcl) const
{
  double res = 0;
  CellType temp = 0;
  size_t valid_pts = 0;
  for (size_t i = 0; i < pcl.size(); ++i) {
    temp = getPtScore(Eigen::Vector2f(pcl.points[i].x, pcl.points[i].y));
    if (temp == OCCUPIED_CELL)
      ++valid_pts;
    res += temp;
  }
  if (pcl.size() == 0)
    return 0;
  return res / (OCCUPIED_CELL * target_points_);
}

template <typename PointType>
double LookUpTable<PointType>::getScore(const std::vector<IndexPoint> &pcl) const
{
  double res = 0;
  CellType temp = 0;
  size_t valid_pts = 0;
  for (size_t i = 0; i < pcl.size(); ++i) {
    if (!(pcl[i].x_idx_ < 0 || pcl[i].y_idx_ < 0 ||
          pcl[i].x_idx_ >= cell_count_row_ ||
          pcl[i].y_idx_ >= cell_count_col_)) {
      ++valid_pts;
    }
    res += getPtScore(pcl[i]);
    // if(temp == OCCUPIED_CELL)
    //  ++valid_pts;
    // res += temp;
  }
  if (pcl.size() == 0)
    return 0;
  // ROS_DEBUG_STREAM("valid points: "<<valid_pts <<" score: "<<res <<"
  // (pts*score)/100 "<<res*valid_pts / 100);
  return res / (target_points_ * OCCUPIED_CELL);
}
template <typename PointType>
double LookUpTable<PointType>::getMaxScore() const
{
  double sum = 0;
  size_t pts = 0;
  for (auto &&cell : table_) {
    sum += cell;
    if (cell > 0)
      ++pts;
  }
  return sum;  /// (pts * OCCUPIED_CELL);
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
void LookUpTable<PointType>::transformIndexes(
    const std::vector<IndexPoint> &source, std::vector<IndexPoint> &out,
    float dx, float dy) const
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

//////////PROTECTED
template <typename PointType>
void LookUpTable<PointType>::initGridWindowFce(
    const pcl::PointCloud<PointType> &target)
{
  table_.clear();
  table_.reserve(cell_count_row_ * cell_count_col_);
  for (size_t i = 0; i < cell_count_row_ * cell_count_col_; ++i) {
    table_.push_back(0);
  }

  // initialize target grid
  // every point is projected to certain cell. If is cell in the grid we set
  // this cell OCCUPIED
  for (size_t i = 0; i < target.size(); ++i) {
    if (target[i].x < maxx_ && target[i].y < maxy_ && target[i].x > minx_ &&
        target[i].y > miny_) {
      Eigen::Vector2f pt(target[i].x, target[i].y);
      size_t idx = getCellIdx(pt);
      table_[idx] = OCCUPIED_CELL;
    }
  }

  // apply smoothing kernel
  for (size_t row = 1; row < cell_count_col_ - 1; ++row) {
    for (size_t col = 1; col < cell_count_row_ - 1; ++col) {
      size_t idx = getCellIdx(row, col);
      if (table_[idx] == OCCUPIED_CELL)
        applyKernel(row, col);
    }
  }
}

template <typename PointType>
size_t LookUpTable<PointType>::getCellIdx(const Eigen::Vector2f &pt) const
{
  size_t x_id = static_cast<size_t>(
      std::floor((pt(0) - minx_ + cell_size_half_) / cell_size_) +
      border_size_);
  size_t y_id = static_cast<size_t>(
      std::floor((-pt(1) + maxy_ + cell_size_half_) / cell_size_) +
      border_size_);
  return y_id * cell_count_row_ + x_id;
}

template <typename PointType>
size_t LookUpTable<PointType>::getCellIdx(size_t row, size_t col) const
{
  return row * cell_count_row_ + col;
}

template <typename PointType>
void LookUpTable<PointType>::applyKernel(size_t row, size_t col)
{
  int half_kernel = kernel_.halfSize();
  for (int r = -half_kernel; r <= half_kernel; ++r) {
    for (int c = -half_kernel; c <= half_kernel; ++c) {
      size_t idx = getCellIdx(row + r, col + c);
      CellType kernel_val = kernel_(r, c);
      if (table_[idx] < kernel_val) {
        table_[idx] = kernel_val;
      }
    }
  }
}

template <typename PointType>
float LookUpTable<PointType>::getPtScore(const Eigen::Vector2f &pt) const
{
  if (pt(0) < maxx_ && pt(1) < maxy_ && pt(0) > minx_ && pt(1) > miny_) {
    return table_[getCellIdx(pt)];
  }
  return 0;
}

template <typename PointType>
float LookUpTable<PointType>::getPtScore(const IndexPoint &pt) const
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
    out << std::setw(5) << o.table_[i];
  }
  return out;
}

///////////////////////////////////// PUBLIC UTILS///////////////////

template <typename PointType>
void rotatePointCloud(const pcl::PointCloud<PointType> &source,
                      pcl::PointCloud<PointType> &outp, float angle)
{
  // x' = x cos f - y sin f
  // y' = y cos f + x sin f
  float si = std::sin(angle);
  float co = std::cos(angle);
  Eigen::MatrixXf rotation(4, 4);
  rotation << co, -si, 0, 0, si, co, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  Eigen::MatrixXf points(4, source.size());
  outp.clear();
  outp.reserve(source.size());
  for (size_t i = 0; i < source.size(); ++i) {
    points.block<4, 1>(0, i) << source[i].x, source[i].y, 0, 0;
    //  PointType pt = source[i];
    // pt.x = source[i].x * co - source[i].y * si;
    // pt.y = source[i].x * si + source[i].y * co;
    // outp.push_back(pt);
  }

  points = rotation * points;
  for (size_t i = 0; i < source.size(); ++i) {
    outp.push_back(PointType(points(0, i), points(1, i), 0));
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
}  // end of pcl namespace

#endif
