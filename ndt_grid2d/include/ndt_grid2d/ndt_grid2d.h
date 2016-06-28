#ifndef NDT_GRID2D_NDT_GRID2D
#define NDT_GRID2D_NDT_GRID2D

#include <ndt_grid2d/voxel_grid2d.h>
#include <ndt_grid2d/types.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/transforms.h>
#include <dynamic_slam_utils/eigen_tools.h>
#include <dynamic_slam_utils/point_cloud_tools.h>

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
  typedef VoxelGrid2D<PointType>::CellPtrVector CellPtrVector;
  typedef VoxelGrid2D<PointType>::Iterator Iterator;
  typedef VoxelGrid2D<PointType>::ConstIterator ConstIterator;

public:
  explicit NDTGrid2D(double timestamp = 0.0);
  // pcl should be in world coordinate frame
  void mergeIn(const PointCloud &pcl, const Pose &origin, bool resize = true);
  void mergeIn(const SelfType &grid, bool transform = true, bool resize = true);
  void mergeIn(const std::vector<CellType> &cells, bool resize);
  void mergeIn(std::vector<CellType> &&cells, bool resize);

  void mergeInTraced(const PointCloud &pcl, const Pose &origin,
                     bool resize = true);
  void mergeInTraced(const SelfType &grid, bool transform, bool resize = true);

  void initialize(const PointCloud &pcl);
  void addPoint(const PointType &pt);
  void computeNDTCells();
  // enlarges grid in all directions to size based on parameters. If parameters
  // are smaller than current state no enlarging or resizing is done.
  void enlarge(float left, float down, float right, float up)
  {
    grid_.enlarge(left, down, right, up);
  }

  // returns only cells with gaussian inside. Includes all cells in radius plus
  // cell where pt belongs.
  CellPtrVector getNeighbors(const Eigen::Vector2d &pt, size_t radius) const;
  bool cellExists(const Point &pt);
  CellType &operator[](const PointType &pt);
  bool isInside(const PointType &pt);
  OccupancyGrid createOccupancyGrid() const;
  // multiple 2 means 2x2 grids will be merged to one cell. Multiple 3  means
  // 4x4 grids will be merged together
  // multiple - multiple of cell_sizes
  SelfType createCoarserGrid(size_t multiple) const;

  void transform(const Eigen::Matrix3d &transform);
  void move(const Eigen::Vector2d &translation);
  void setCellSize(float cell_size)
  {
    cell_size_ = cell_size;
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

protected:
  Eigen::Vector3d origin_;  // x,y,theta in the world frame
  float cell_size_;
  bool initialized_;
  double timestamp_;
  DataGrid grid_;

  SelfType createGrid(const PointCloud &pcl) const;
  void transformNDTCells(std::vector<CellType> &grid,
                         const Eigen::Matrix3d &transform);
};

// //////////////////IMPLEMENTATION ///////////
template <typename CellType, typename PointType>
NDTGrid2D<CellType, PointType>::NDTGrid2D(double timestamp)
  : origin_(Eigen::Vector3d::Zeros())
  , cell_size_(0.25)
  , initialized_(false)
  , timestamp_(timestamp)
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
    transformNDTCells(occupied_cells, trans);
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
    initialized = true;
    float minx, miny, maxx, maxy;
    getMinMaxNDT2D(cells, minx, miny, maxx, maxy);
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
    initialized = true;
    float minx, miny, maxx, maxy;
    getMinMaxNDT2D(cells, minx, miny, maxx, maxy);
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
    transformNDTCells(occupied_cells, trans);
  }
  CellPtrVector trace_cells;
  for (auto &&cell : occupied_cells) {
    if (cell.hasGaussian()) {
      CellPtrVector ray_cells = grid_.rayTrace(cell.getMean());
      for (CellType *ray_cell : ray_cells) {
        // update cell on line between start and end point
        ray_cell->updateOccupancy(this->origin_, cell.getMean(), cell.points());
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
void NDTGrid2D<CellType, PointType>::addPoint(const PointType &pt)
{
  DataGrid::Point point;
  point << pt.x, pt.y;
  CellType *cell = grid_.getCellPtr(point);
  if (!existsCell(point)) {
    // this field doesn't exist in cell yet
    grid_.addCell(point, CellType(), false);
    grid_.getCellPtr(point)->addPoint(point);
  } else {
    grid_[point]->addPoint(point);
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
NDTGrid2D<CellType, PointType>::CellPtrVector
NDTGrid2D<CellType, PointType>::getNeighbors(const PointType &pt,
                                             size_t radius) const
{
  GridData::Point point;
  point << pt.x, pt.y;
  CellPtrVector neighbours = grid_.getNeighbors(point, radius);
  CellPtrVector gaussians;
  gaussians.resize(radius + 1 * 4);
  for (CellType *cell : neighbours) {
    if (cell->hasGaussian())
      gaussians.push_back(cell);
  }
  return gaussians;
}

template <typename CellType, typename PointType>
bool NDTGrid2D<CellType, PointType>::cellExists(const Point &pt)
{
  GridData::Point point;
  point << pt.x, pt.y;
  return grid_.cellExists(point);
}

template <typename CellType, typename PointType>
CellType &NDTGrid2D<CellType, PointType>::operator[](const PointType &pt)
{
  GridData::Point point;
  point << pt.x, pt.y;
  return grid_[point];
}

template <typename CellType, typename PointType>
bool NDTGrid2D<CellType, PointType>::isInside(const PointType &pt)
{
  GridData::Point point;
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
  GridData::CellPtrVector valid_cells = grid_.getValidCellsPtr();
  for (auto cell_it : grid_) {
    if (*cell_it == nullptr)
      oc_grid.cells_.push_back(-1);
    else
      oc_grid.cells_.push_back(cell_it->getOccupancy());
  }
  return oc_grid;
}

template <typename CellType, typename PointType>
NDTGrid2D<CellType, PointType>::SelfType
NDTGrid2D<CellType, PointType>::createCoarserGrid(size_t multiple) const
{
  if (multiple <= 1)
    return grid_;

  SelfType coarse_grid;
  coarse_grid.setCellSize(multiple * cell_size_);
  mergeIn(grid_.getValidCells(), true);
  return coarse_grid;
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::transform(const Eigen::Matrix3d &transform)
{
  // get valid cells apply on them transformation and fuse them back to empty
  // grid
  DataGrid trans_grid = grid_.clone();
  CellVector cells = grid_.getValidCells();
  transformNDTCells(cells, transform);
  trans_grid.mergeIn(cells, true);
  grid_ = std::move(trans_grid);
  // transform origin of ndt grid
  eigt::transform2d_t trans_mat(transform);
  origin_ = trans_mat * origin_;
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::move(const Eigen::Vector2d &translation)
{
  // TODO TODO TODO change calculation. Take in respect real movement
  // grid_.translate(translation,true);
  // origin_(0) += translation(0)* cell_size_;
  // origin_(1) += translation(1) * cell_size_;
}

////////////////////PROTECCTED////////////////
template <typename CellType, typename PointType>
NDTGrid2D<CellType, PointType>::SelfType
NDTGrid2D<CellType, PointType>::createGrid(const PointCloud &pcl) const
{
  SelfType localGrid;
  localGrid.setCellSize(cell_size_);
  localGrid.setOrigin(this->origin_);
  float minx, miny, maxx, maxy;
  pcl::getMinMax2D(trans_pcl, minx, miny, maxx, maxy);
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
