#ifndef GRAPH_SLAM_UK_NDT_GRID2D
#define GRAPH_SLAM_UK_NDT_GRID2D

#include <graph_slam_uk/ndt/voxel_grid2d.h>
#include <graph_slam_uk/ndt/ndt_grid2d_interface.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/transforms.h>
#include <graph_slam_uk/utils/eigen_tools.h>
#include <graph_slam_uk/utils/point_cloud_tools.h>
#include <graph_slam_uk/registration/scanmatchable_grid_interface.h>
#include <boost/shared_ptr.hpp>
#include <graph_slam_uk/ndt/output_msgs.h>

namespace slamuk
{
template <typename CellType, typename PointType>
class NDTGrid2D
{
protected:
  typedef VoxelGrid2D<CellType> DataGrid;

public:
  typedef pcl::PointCloud<PointType> PointCloud;
  typedef NDTGrid2D<CellType, PointType> SelfType;
  typedef CellType Cell;
  typedef Eigen::Vector3d Pose;
  typedef typename VoxelGrid2D<CellType>::CellPtrVector CellPtrVector;
  typedef typename VoxelGrid2D<CellType>::Iterator Iterator;
  typedef typename VoxelGrid2D<CellType>::ConstIterator ConstIterator;
  typedef boost::shared_ptr<SelfType> Ptr;
  typedef boost::shared_ptr<const SelfType> ConstPtr;

public:
  explicit NDTGrid2D(double timestamp = 0.0);
  /////////////IMPLEMENTATION oF INDTGrid2D
  // pcl should be in world coordinate frame
  void initialize(const PointCloud &pcl);
  void initializeSimple(const PointCloud &pcl);
  void mergeIn(const PointCloud &pcl, const Pose &origin, bool resize = true);
  void mergeIn(const SelfType &grid, bool transform = true, bool resize = true);
  void mergeIn(const std::vector<CellType> &cells, bool resize);

  void mergeInTraced(const PointCloud &pcl, const Pose &origin,
                     bool resize = true);
  void mergeInTraced(const SelfType &grid, bool transform, bool resize = true);

  // enlarges grid in all directions to size based on parameters. If parameters
  // are smaller than current state no enlarging or resizing is done.
  void enlarge(float left, float down, float right, float up)
  {
    grid_.enlarge(left, down, right, up);
  }
  void move(const Eigen::Vector2d &translation);
  void transform(const Eigen::Matrix3d &transform);
  OccupancyGrid createOccupancyGrid() const;

  void setCellSize(float cell_size)
  {
    cell_size_ = cell_size;
  }
  float getCellSize() const
  {
    return cell_size_;
  }

  Eigen::Vector3d getOrigin() const
  {
    return origin_;
  }
  void setOrigin(const Eigen::Vector3d &origin)
  {
    origin_ = origin;
  }

  double getTimestamp() const
  {
    return timestamp_;
  }
  void setTimestamp(double timestamp)
  {
    timestamp_ = timestamp;
  }
  bool operator<(const NDTGrid2D &other) const
  {
    return timestamp_ < other.timestamp_;
  }
  /////////////// IMPLEMENTATION of IScanmatchableGrid
  // returns only cells with gaussian inside. Includes all cells in radius plus
  // cell where pt belongs.
  CellPtrVector getNeighbors(const Eigen::Vector2d &pt, size_t radius) const;

  // multiple 2 means 2x2 grids will be merged to one cell. Multiple 3  means
  // 4x4 grids will be merged together
  // multiple - multiple of cell_sizes
  SelfType createCoarserGrid(size_t multiple) const;

  PointCloud getMeans() const;

  CellPtrVector getGaussianCells() const;

  /////////////// ADITIONAL
  CellType &operator[](const PointType &pt);
  bool isInside(const PointType &pt);
  Ptr makeShared() const
  {
    return Ptr(new SelfType(*this));
  }

  double getRadius() const
  {
    return std::max(grid_.width() * cell_size_, grid_.height() * cell_size_);
  }
  Eigen::Vector2d getCentroid() const
  {
    return origin_.head(2);
  }
  NDTGridMsg serialize() const;

protected:
  Eigen::Vector3d origin_;  // x,y,theta in the world frame
  float cell_size_;
  bool initialized_;
  double timestamp_;
  DataGrid grid_;
  Eigen::Vector2d cumul_move_;

  void mergeIn(std::vector<CellType> &&cells, bool resize);
  void addPoint(const PointType &pt);
  void computeNDTCells();
  SelfType createGrid(const PointCloud &pcl) const;
  void transformNDTCells(std::vector<CellType> &grid,
                         const Eigen::Matrix3d &transform);
};

// //////////////////IMPLEMENTATION ///////////
template <typename CellType, typename PointType>
NDTGrid2D<CellType, PointType>::NDTGrid2D(double timestamp)
  : origin_(Eigen::Vector3d::Zero())
  , cell_size_(0.25)
  , initialized_(false)
  , timestamp_(timestamp)
  , grid_()
  , cumul_move_(0, 0)
{
  grid_.setCellSize(cell_size_);
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::mergeIn(const PointCloud &pcl,
                                             const Pose &origin, bool resize)
{
  eigt::transform2d_t<double> trans =
      eigt::transBtwPoses(this->origin_, origin);
  PointCloud trans_pcl;
  // transforming input pcl to reference frame of this grid
  pcl::transformPointCloud(pcl, trans_pcl,
                           eigt::convertFromTransform(trans).cast<float>());

  mergeIn(createGrid(trans_pcl), false, resize);
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::mergeIn(const SelfType &grid,
                                             bool transform, bool resize)
{
  // all valid cells. Including cells without proper gaussian
  std::vector<CellType> occupied_cells = grid.grid_.getValidCells();
  if (transform) {
    // we need to transform each cell to its new position
    eigt::transform2d_t<double> trans =
        eigt::transBtwPoses(this->origin_, grid.origin_);
    transformNDTCells(occupied_cells, trans.matrix());
  }
  mergeIn(occupied_cells, resize);
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::mergeIn(const std::vector<CellType> &cells,
                                             bool resize)
{
  // incoming vector of cells includes all cells [cells with gaussian, visited
  // unoccupied cells,unoccupied cells with gaussian]
  // only cells with gaussian are used
  if (resize || !initialized_) {
    initialized_ = true;
    float minx, miny, maxx, maxy;
    pcl::getMinMaxNDT2D(cells, minx, miny, maxx, maxy);
    grid_.enlarge(minx, miny, maxx, maxy);
    for (auto &&cell : cells) {
      if (cell.hasGaussian())
        grid_.addCell(cell.getMean(), cell);
    }
  } else {
    for (auto &&cell : cells) {
      // cell is discarted if it is outside of the boundries for current grid_
      if (cell.hasGaussian())
        if (grid_.isInside(cell.getMean()))
          grid_.addCell(cell.getMean(), cell);
    }
  }
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::mergeIn(std::vector<CellType> &&cells,
                                             bool resize)
{
  // incoming vector of cells includes all cells [cells with gaussian, visited
  // unoccupied cells,unoccupied cells with gaussian]
  if (resize || !initialized_) {
    initialized_ = true;
    float minx, miny, maxx, maxy;
    pcl::getMinMaxNDT2D(cells, minx, miny, maxx, maxy);
    grid_.enlarge(minx, miny, maxx, maxy);
    for (auto &&cell : cells) {
      grid_.addCell(cell.getMean(), std::move(cell));
    }
  } else {
    for (auto &&cell : cells) {
      // cell is discarted if it is outside of the boundries for current grid_
      if (grid_.isInside(cell.getMean()))
        grid_.addCell(cell.getMean(), std::move(cell));
    }
  }
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::mergeInTraced(const PointCloud &pcl,
                                                   const Pose &origin,
                                                   bool resize)
{
  eigt::transform2d_t<double> trans =
      eigt::transBtwPoses(this->origin_, origin);
  PointCloud trans_pcl;
  // transforming input pcl to reference frame of this grid
  pcl::transformPointCloud(pcl, trans_pcl,
                           eigt::convertFromTransform(trans).cast<float>());
  mergeInTraced(createGrid(trans_pcl), false, resize);
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::mergeInTraced(const SelfType &grid,
                                                   bool transform, bool resize)
{
  // all valid cells. Including cells without proper gaussian
  std::vector<CellType> occupied_cells = grid.grid_.getValidCells();
  if (transform) {
    // we need to transform each cell to its new position
    eigt::transform2d_t<double> trans =
        eigt::transBtwPoses(this->origin_, grid.origin_);
    transformNDTCells(occupied_cells, trans.matrix());
  }
  for (auto &&cell : occupied_cells) {
    if (cell.hasGaussian()) {
      CellPtrVector ray_cells = grid_.rayTrace(cell.getMean());
      for (CellType *ray_cell : ray_cells) {
        // update cell on line between start and end point
        ray_cell->updateOccupancy(this->origin_.head(2), cell.getMean(),
                                  cell.points());
      }
    }
  }
  mergeIn(grid.grid_.getValidCells(), resize);
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::initialize(const PointCloud &pcl)
{
  grid_.clear();
  mergeInTraced(pcl, origin_, true);
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::initializeSimple(const PointCloud &pcl)
{
  grid_.clear();
  mergeIn(pcl, origin_, true);
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::addPoint(const PointType &pt)
{
  typename DataGrid::Point point;
  point << pt.x, pt.y;
  if (!grid_.cellExists(point)) {
    // this field doesn't exist in cell yet
    grid_.addCell(point, CellType(), false);
    grid_[point].addPoint(point);
  } else {
    grid_[point].addPoint(point);
  }
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::computeNDTCells()
{
  for (CellType *cell : grid_.getValidCellsPtr()) {
    cell->computeGaussian();
  }
}

template <typename CellType, typename PointType>
typename NDTGrid2D<CellType, PointType>::CellPtrVector
NDTGrid2D<CellType, PointType>::getNeighbors(const Eigen::Vector2d &pt,
                                             size_t radius) const
{
  CellPtrVector neighbours = grid_.getNeighbors(pt, radius);
  CellPtrVector gaussians;
  gaussians.reserve(neighbours.size());
  for (CellType *cell : neighbours) {
    if (cell->hasGaussian())
      gaussians.push_back(cell);
  }
  return gaussians;
}

template <typename CellType, typename PointType>
CellType &NDTGrid2D<CellType, PointType>::operator[](const PointType &pt)
{
  typename DataGrid::Point point;
  point << pt.x, pt.y;
  return grid_[point];
}

template <typename CellType, typename PointType>
bool NDTGrid2D<CellType, PointType>::isInside(const PointType &pt)
{
  typename DataGrid::Point point;
  point << pt.x, pt.y;
  return grid_.isInside(point);
}

template <typename CellType, typename PointType>
OccupancyGrid NDTGrid2D<CellType, PointType>::createOccupancyGrid() const
{
  OccupancyGrid oc_grid;
  oc_grid.origin_ = origin_;
  oc_grid.width_ = grid_.width();
  oc_grid.height_ = grid_.height();
  oc_grid.resolution_ = cell_size_;
  oc_grid.cells_.resize(oc_grid.width_ * oc_grid.height_);
  typename DataGrid::CellPtrVector valid_cells = grid_.getValidCellsPtr();
  for (CellType *cell_it : valid_cells) {
    if (cell_it == nullptr)
      oc_grid.cells_.push_back(-1);
    else
      oc_grid.cells_.push_back(cell_it->getOccupancy());
  }
  return oc_grid;
}

template <typename CellType, typename PointType>
typename NDTGrid2D<CellType, PointType>::SelfType
NDTGrid2D<CellType, PointType>::createCoarserGrid(size_t multiple) const
{
  if (multiple <= 1)
    return *this;

  SelfType coarse_grid;
  coarse_grid.setCellSize(multiple * cell_size_);
  coarse_grid.mergeIn(grid_.getValidCells(), true);
  return coarse_grid;
}

template <typename CellType, typename PointType>
typename NDTGrid2D<CellType, PointType>::PointCloud
NDTGrid2D<CellType, PointType>::getMeans() const
{
  PointCloud pcl;
  for (auto &&cell : getGaussianCells()) {
    PointType pt;
    Eigen::Vector2d mean = cell->getMean();
    pt.x = mean(0);
    pt.y = mean(1);
    pt.z = 0;
  }
}

template <typename CellType, typename PointType>
typename NDTGrid2D<CellType, PointType>::CellPtrVector
NDTGrid2D<CellType, PointType>::getGaussianCells() const
{
  CellPtrVector res;
  res.reserve(grid_.validCells());
  for (auto &&cell : grid_.getValidCellsPtr()) {
    if (cell->hasGaussian())
      res.push_back(cell);
  }
  return res;
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::transform(const Eigen::Matrix3d &transform)
{
  // get valid cells apply on them transformation and fuse them back to empty
  // grid
  DataGrid trans_grid = grid_.clone();
  typename DataGrid::CellVector cells = grid_.getValidCells();
  transformNDTCells(cells, transform);
  trans_grid.mergeIn(cells, true);
  grid_ = std::move(trans_grid);
  // transform origin of ndt grid
  eigt::transform2d_t<double> trans_mat(transform);
  origin_ = trans_mat * origin_;
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::move(const Eigen::Vector2d &translation)
{
  cumul_move_ += translation;
  Eigen::Vector2i move(0, 0);
  if (std::abs(std::floor(cumul_move_(0) / cell_size_)) > 0) {
    move(0) = std::floor(cumul_move_(0) / cell_size_);
    cumul_move_(0) -= move(0);
  }
  if (std::abs(std::floor(cumul_move_(1) / cell_size_)) > 0) {
    move(1) = std::floor(cumul_move_(1) / cell_size_);
    cumul_move_(1) -= move(1);
  }
  // move underling voxel grid
  grid_.translate(move, true);
  // move origin of this grid
  origin_(0) += move(0) * cell_size_;
  origin_(1) += move(1) * cell_size_;
  // applies transformation to all means of ndt cells
  for (auto &&cell : getGaussianCells()) {
    cell->getMean() += move.cast<double>();
  }
}

template <typename CellType, typename PointType>
NDTGridMsg NDTGrid2D<CellType, PointType>::serialize() const
{
  NDTGridMsg msg;
  msg.origin_ << origin_(0), origin_(1), 0;
  msg.size_ << grid_.width() * cell_size_, grid_.height() * cell_size_, 0;
  msg.cell_sizes_ << cell_size_, cell_size_, 0;
  for (auto &&cell : grid_.getValidCellsPtr) {
    msg.cells_.push_back(cell->serialize());
  }
  return msg;
}

////////////////////PROTECCTED////////////////
template <typename CellType, typename PointType>
typename NDTGrid2D<CellType, PointType>::SelfType
NDTGrid2D<CellType, PointType>::createGrid(const PointCloud &pcl) const
{
  SelfType localGrid;
  localGrid.setCellSize(cell_size_);
  localGrid.setOrigin(this->origin_);
  float minx, miny, maxx, maxy;
  pcl::getMinMax2D(pcl, minx, miny, maxx, maxy);
  localGrid.enlarge(minx, miny, maxx, maxy);
  for (size_t i = 0; i < pcl.size(); ++i) {
    if (std::isnan(pcl[i].x) || std::isnan(pcl[i].y) || std::isnan(pcl[i].z))
      continue;
    localGrid.addPoint(pcl[i]);
  }
  localGrid.computeNDTCells();
  return localGrid;
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::transformNDTCells(
    std::vector<CellType> &grid, const Eigen::Matrix3d &transform)
{
  for (CellType &cell : grid) {
    cell.transform(transform);
  }
}

}  // end of namespace slamuk

#endif
