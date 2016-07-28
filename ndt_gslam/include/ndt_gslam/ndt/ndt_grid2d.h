#ifndef NDT_GSLAM_NDT_GRID2D
#define NDT_GSLAM_NDT_GRID2D

#include <ndt_gslam/ndt/output_msgs.h>
#include <ndt_gslam/ndt/voxel_grid2d.h>
#include <ndt_gslam/utils/eigen_tools.h>
#include <ndt_gslam/utils/point_cloud_tools.h>
#include <pcl/common/transforms.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <boost/shared_ptr.hpp>

namespace slamuk
{
/**
 * @brief      Datastructure representing NDT grid.
 *
 * @tparam     CellType   Type of NDT cell
 * @tparam     PointType  Type of PCL point type used for initialization.
 */
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
  typedef Eigen::Transform<double, 2, Eigen::TransformTraits::Affine> Transform;
  typedef typename VoxelGrid2D<CellType>::CellPtrVector CellPtrVector;
  typedef typename VoxelGrid2D<CellType>::Iterator Iterator;
  typedef typename VoxelGrid2D<CellType>::ConstIterator ConstIterator;
  typedef boost::shared_ptr<SelfType> Ptr;
  typedef boost::shared_ptr<const SelfType> ConstPtr;

public:
  explicit NDTGrid2D(double timestamp = 0.0);
  NDTGrid2D(const Eigen::Vector3d &origin, double timestamp = 0.0);

  /**
   * @brief      Creates NDT grid from point cloud and update occupancy of empty
   *             spaces.
   *
   * @param[in]  pcl   The pcl
   */
  void initialize(const PointCloud &pcl);

  /**
   * @brief      Creates NDT grid from point cloud without occupancy update.
   *
   * @param[in]  pcl   The pcl
   */
  void initializeSimple(const PointCloud &pcl);

  /**
   * @brief      Merges in point cloud based on it origin in wolrd coordiate
   *             system.
   *
   * @param[in]  pcl     The pcl
   * @param[in]  origin  The origin of the point cloud
   * @param[in]  resize  Set ture if you want to enlarge grid based to fit
   *                     inputed point cloud
   */
  void mergeIn(const PointCloud &pcl, const Pose &origin, bool resize = true);

  /**
   * @brief      Merge in another NDt grid.
   *
   * @param[in]  grid       The grid
   * @param[in]  transform  Set true if you want to apply change bases
   *                        transformation to all cells in inserted grid
   * @param[in]  resize     Set ture if you want to enlarge grid based to fit
   *                     inputed NDT grid's cells
   */
  void mergeIn(const SelfType &grid, bool transform = true, bool resize = true);

  /**
   * @brief      Merge in all cells from vector.
   *
   * @param[in]  cells   The cells
   * @param[in]  resize  Set ture if you want to enlarge grid base to fit
   *                     all cells.
   */
  void mergeIn(const std::vector<CellType> &cells, bool resize);

  /**
   * @brief      Merges in point cloud based on it origin in wolrd coordiate
   *             system. Apply occupancy update to identifie empty space.
   *
   * @param[in]  pcl     The pcl
   * @param[in]  origin  The origin of the point cloud
   * @param[in]  resize  Set ture if you want to enlarge grid based to fit
   *                     inputed point cloud
   */
  void mergeInTraced(const PointCloud &pcl, const Pose &origin,
                     bool resize = true);

  /**
   * @brief      Merge in another NDt grid. Apply occupancy update to identifie
   * empty space.
   *
   * @param[in]  grid       The grid
   * @param[in]  transform  Set true if you want to apply change bases
   *                        transformation to all cells in inserted grid
   * @param[in]  resize     Set ture if you want to enlarge grid based to fit
   *                     inputed NDT grid's cells
   */
  void mergeInTraced(const SelfType &grid, bool transform, bool resize = true);

  // enlarges grid in all directions to size based on parameters. If parameters
  // are smaller than current state no enlarging or resizing is done.
  //
  // @param[in]  left   The left
  // @param[in]  down   The down
  // @param[in]  right  The right
  // @param[in]  up     { parameter_description }
  //

  /**
   * @brief      Enlarges grid in all directions to size based on parameters. If
   *             parameters are smaller than current state no enlarging or
   *             resizing is done.
   *
   * @param[in]  left   The left direction
   * @param[in]  down   The down direction
   * @param[in]  right  The right direction
   * @param[in]  up     The up direction
   */
  void enlarge(float left, float down, float right, float up)
  {
    grid_.enlarge(left, down, right, up);
  }

  /**
   * @brief      moves grid in horizontal or vertical direction in units of
   *             whole cell. Transformation which has not been used is returnet.
   *             Always keeps this grid alligned whith Kartesian coordinates and
   *             doesn't change cells covariances. It is faster and more robust
   *             than regular transform method.
   *
   * @param[in]  transform  transformation which should be applied to this grid
   *
   * @return     residual (unused) transformation
   */
  Transform move(const Transform &transform);

  /**
   * @brief      Transforms whole grid.
   *
   * @param[in]  transform  The transform
   */
  void transform(const Transform &transform);

  /**
   * @brief      Creates an occupancy grid message.
   *
   * @return     The occupancy message.
   */
  OccupancyGrid createOccupancyGrid() const;

  /**
   * @brief      Sets the cell size.
   *
   * @param[in]  cell_size  The cell size
   */
  void setCellSize(float cell_size)
  {
    cell_size_ = cell_size;
    grid_.setCellSize(cell_size);
  }

  /**
   * @brief      Gets the cell size.
   *
   * @return     The cell size.
   */
  float getCellSize() const
  {
    return cell_size_;
  }

  /**
   * @brief      Gets the origin.
   *
   * @return     The origin.
   */
  Eigen::Vector3d getOrigin() const
  {
    return origin_;
  }

  /**
   * @brief      Sets the origin.
   *
   * @param[in]  origin  The origin
   */
  void setOrigin(const Eigen::Vector3d &origin)
  {
    origin_ = origin;
  }

  /**
   * @brief      Gets the time stamp.
   *
   * @return     The timestamp.
   */
  double getTimestamp() const
  {
    return timestamp_;
  }
  /**
   * @brief      Sets the time-stamp.
   *
   * @param[in]  timestamp  The time-stamp
   */
  void setTimestamp(double timestamp)
  {
    timestamp_ = timestamp;
  }

  /**
   * @brief      The comparator based on time-stamp of creation.
   *
   * @param[in]  other  The other grid
   *
   * @return     true if current grid is smaller than other grid
   */
  bool operator<(const NDTGrid2D &other) const
  {
    return timestamp_ < other.timestamp_;
  }

  // returns only cells with Gaussian inside. Includes all cells in radius plus
  // cell where pt belongs.

  /**
   * @brief      The radius search around certain point. Returns only cells with
   *             Gaussian inside. Includes all cells in radius plus cell where
   *             pt belongs.
   *
   * @param[in]  pt      The point
   * @param[in]  radius  The radius
   *
   * @return     The neighbors.
   */
  CellPtrVector getNeighbors(const Eigen::Vector2d &pt, float radius) const;

  /**
   * @brief      Gets the k nearest neighbors with Gaussian information inside.
   *
   * @param[in]  pt    The point
   * @param[in]  K     The number of nearest neighbors
   *
   * @return     The k nearest neighbors.
   */
  CellPtrVector getKNearestNeighbors(const Eigen::Vector2d &pt, size_t K) const;

  /**
   * @brief      Creates a coarser grid based on selected cell size. Apply
   *             merging of cells. It is advised to use multiples of 2 from the
   *             base cell size.
   *
   * @param[in]  cell_size  The cell size
   *
   * @return     { description_of_the_return_value }
   */
  SelfType createCoarserGrid(float cell_size) const;

  /**
   * @brief      Gets the point cloud with all means in the grid
   *
   * @return     The means.
   */
  const typename PointCloud::Ptr getMeans() const
  {
    return means_;
  }

  /**
   * @brief      Gets the means transformed to from local grid centric
   * coordinate
   *             system into global coordinate system.
   *
   * @return     The means in global frame.
   */
  const typename PointCloud::Ptr getMeansTransformed() const
  {
    typename PointCloud::Ptr res(new PointCloud());
    Eigen::Matrix4d trans = eigt::convertFromTransform(
        eigt::transBtwFrames(origin_, Eigen::Vector3d(0, 0, 0)));
    pcl::transformPointCloud(*means_, *res, trans);
    return res;
  }

  /**
   * @brief      Gets all cells containing Gaussian distribution.
   *
   * @return     The gaussian cells.
   */
  CellPtrVector getGaussianCells() const;

  /////////////// ADITIONAL

  /**
   * @brief      Returns cell occupying cell at point.
   *
   * @param[in]  pt    The point
   *
   * @return     The cell.
   */
  const CellType &operator[](const PointType &pt) const;

  /**
     * @brief      Returns cell occupying cell at point.
     *
     * @param[in]  pt    The point
     *
     * @return     The cell.
     */
  const CellType &operator[](const Eigen::Vector2d &pt) const;

  /**
     * @brief      Returns cell occupying cell at point.
     *
     * @param[in]  pt    The point
     *
     * @return     The cell.
     */
  CellType &operator[](const PointType &pt);
  /**
   * @brief      Returns cell occupying cell at point.
   *
   * @param[in]  pt    The point
   *
   * @return     The cell.
   */
  CellType &operator[](const Eigen::Vector2d &pt);

  /**
   * @brief      Determines if point is inside the bounds of grid.
   *
   * @param[in]  pt    The point
   *
   * @return     True if inside, False otherwise.
   */
  bool isInside(const PointType &pt);
  /**
 * @brief      Determines if point is inside the bounds of grid.
 *
 * @param[in]  pt    The point
 *
 * @return     True if inside, False otherwise.
 */
  bool isInside(const Eigen::Vector2d &pt) const;

  /**
   * @brief      Makes a shared pointer from this grid.
   *
   * @return     The shared pointer.
   */
  Ptr makeShared() const
  {
    return Ptr(new SelfType(*this));
  }

  /**
   * @brief      Gets the radius of farthest cells from centroid of this grid.
   *
   * @return     The radius.
   */
  double getRadius() const
  {
    return grid_.width() * cell_size_ + grid_.height() * cell_size_ / 4;
  }
  /**
   * @brief      Gets the centroid of the grid.
   *
   * @return     The centroid.
   */
  Eigen::Vector2d getCentroid() const
  {
    double x = grid_.right() * cell_size_ - grid_.left() * cell_size_;
    double y = (grid_.up() * cell_size_ - grid_.down() * cell_size_);
    return Eigen::Vector2d(x, y);
  }

  /**
   * @brief      Removes all cells. Capacity and size stays the same
   */
  void clear()
  {
    grid_.clear();
    initialized_ = false;
  }

  /**
   * @brief      Grid to message serialization.
   *
   * @return     message
   */
  NDTGridMsg serialize() const;

  template <typename Cell, typename Pt>
  friend std::ostream &operator<<(std::ostream &os,
                                  const NDTGrid2D<Cell, Pt> &grid);

protected:
  Eigen::Vector3d origin_;  // x,y,theta in the world frame
  float cell_size_;
  bool initialized_;
  double timestamp_;
  DataGrid grid_;
  pcl::KdTreeFLANN<PointType> kdtree_;
  typename PointCloud::Ptr means_;

  void mergeIn(std::vector<CellType> &&cells, bool resize);
  void addPoint(const PointType &pt);
  void computeNDTCells();
  SelfType createGrid(const PointCloud &pcl) const;
  void transformNDTCells(std::vector<CellType> &grid,
                         const Transform &transform);
  void updateKDTree()
  {
    means_ = getMeans();
    kdtree_.setInputCloud(means_, nullptr);
  }
  CellPtrVector getKNearestNeighborsVoxel(const Eigen::Vector2d &pt,
                                          size_t K) const;
  CellPtrVector getKNearestNeighborsKDTree(const Eigen::Vector2d &pt,
                                           size_t K) const;
  void updateMeansCloud();
};

// //////////////////IMPLEMENTATION ///////////
template <typename CellType, typename PointType>
NDTGrid2D<CellType, PointType>::NDTGrid2D(double timestamp)
  : origin_(Eigen::Vector3d::Zero())
  , cell_size_(0.25)
  , initialized_(false)
  , timestamp_(timestamp)
  , grid_()
  , kdtree_()
  , means_(new PointCloud())
{
  grid_.setCellSize(cell_size_);
}

template <typename CellType, typename PointType>
NDTGrid2D<CellType, PointType>::NDTGrid2D(const Eigen::Vector3d &origin,
                                          double timestamp)
  : origin_(origin)
  , cell_size_(0.25)
  , initialized_(false)
  , timestamp_(timestamp)
  , grid_()
{
  grid_.setCellSize(cell_size_);
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::mergeIn(const PointCloud &pcl,
                                             const Pose &origin, bool resize)
{
  eigt::transform2d_t<double> trans =
      eigt::transBtwFrames(origin, this->origin_);
  // std::cout << trans.matrix() << std::endl;
  PointCloud trans_pcl;
  // transforming input pcl to reference frame of this grid
  pcl::transformPointCloud(pcl, trans_pcl, eigt::convertFromTransform(trans));
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
        eigt::transBtwFrames(grid.origin_, this->origin_);
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
    initialized_ = true;
    float minx, miny, maxx, maxy;
    pcl::getMinMaxNDT2D(cells, minx, miny, maxx, maxy);
    grid_.enlarge(minx, miny, maxx, maxy);
    for (auto &&cell : cells) {
      if (cell.hasGaussian()) {
        grid_.addCell(cell.getMean().head(2), cell, false);
      }
    }
  } else {
    for (auto &&cell : cells) {
      // cell is discarted if it is outside of the boundries for current grid_
      if (cell.hasGaussian())
        if (grid_.isInside(cell.getMean().head(2)))
          grid_.addCell(cell.getMean().head(2), cell);
    }
  }
  updateMeansCloud();
  updateKDTree();
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::mergeIn(std::vector<CellType> &&cells,
                                             bool resize)
{
  // incoming vector of cells includes all cells [cells with Gaussian, visited
  // unoccupied cells,unoccupied cells with Gaussian]
  if (resize || !initialized_) {
    initialized_ = true;
    float minx, miny, maxx, maxy;
    pcl::getMinMaxNDT2D(cells, minx, miny, maxx, maxy);
    grid_.enlarge(minx, miny, maxx, maxy);
    for (auto &&cell : cells) {
      grid_.addCell(cell.getMean().head(2), std::move(cell));
    }
  } else {
    for (auto &&cell : cells) {
      // cell is discarded if it is outside of the boundaries for current grid_
      if (grid_.isInside(cell.getMean().head(2)))
        grid_.addCell(cell.getMean().head(2), std::move(cell));
    }
  }
  updateMeansCloud();
  updateKDTree();
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::mergeInTraced(const PointCloud &pcl,
                                                   const Pose &origin,
                                                   bool resize)
{
  eigt::transform2d_t<double> trans =
      eigt::transBtwFrames(origin, this->origin_);
  PointCloud trans_pcl;
  // transforming input pcl to reference frame of this grid
  pcl::transformPointCloud(pcl, trans_pcl, eigt::convertFromTransform(trans));
  auto grid = createGrid(trans_pcl);
  grid.setOrigin(origin);
  mergeInTraced(std::move(grid), false, resize);
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::mergeInTraced(const SelfType &grid,
                                                   bool transform, bool resize)
{
  // all valid cells. Including cells without proper gaussian
  std::vector<CellType> occupied_cells = grid.grid_.getValidCells();
  eigt::transform2d_t<double> trans =
      eigt::transBtwFrames(grid.origin_, this->origin_);
  Eigen::Vector2d origin2 = trans * Eigen::Vector2d(0, 0);
  if (transform) {
    // we need to transform each cell to its new position
    transformNDTCells(occupied_cells, trans);
  }
  // conversion from x y theta pose to x,y,z coordinates for 3dCell interface
  typename CellType::Vector start(origin2(0), origin2(1), 0);
  // std::cout << "start: " << start.transpose() << std::endl;
  // go over all cells in grid, whichh is merging in.
  for (auto &&cell : occupied_cells) {
    if (cell.hasGaussian()) {
      // request all cells on ray passing through original grid with start and
      // end of raytracing taken from new merge in grid
      CellPtrVector ray_cells =
          this->grid_.rayTrace(start.head(2), cell.getMean().head(2));
      // std::cout << ray_cells.size() << std::endl;
      for (CellType *ray_cell : ray_cells) {
        // update cell on line between start and end point
        ray_cell->updateOccupancy(start, cell.getMean(), cell.points());
      }
    }
  }
  mergeIn(occupied_cells, resize);
  // std::cout << *this << std::endl << std::endl;
  // char a;
  // std::cin >> a;
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
  Eigen::Vector3d point(pt.x, pt.y, pt.z);
  Eigen::Vector2d point2d = point.head(2);
  if (!grid_.cellExists(point2d)) {
    // this field doesn't exist in cell yet
    grid_.addCell(point2d, CellType(), false);
    grid_[point2d].addPoint(point);
  } else {
    grid_[point2d].addPoint(point);
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
                                             float radius) const
{
  size_t manh_radius;
  CellPtrVector gaussians;
  if (radius <= 0) {
    manh_radius = 0;
    CellPtrVector neighbours = grid_.getNeighbors(pt, manh_radius);
    for (CellType *cell : neighbours) {
      if (cell->hasGaussian())
        gaussians.push_back(cell);
    }
  } else {
    manh_radius = std::ceil(radius / cell_size_);
    CellPtrVector neighbours = grid_.getNeighbors(pt, manh_radius);
    gaussians.reserve(neighbours.size());
    for (CellType *cell : neighbours) {
      if (cell->hasGaussian())
        if ((cell->getMean().head(2) - pt).norm() <= radius)
          gaussians.push_back(cell);
    }
  }
  return gaussians;
}

template <typename CellType, typename PointType>
typename NDTGrid2D<CellType, PointType>::CellPtrVector
NDTGrid2D<CellType, PointType>::getKNearestNeighbors(const Eigen::Vector2d &pt,
                                                     size_t K) const
{
  return getKNearestNeighborsKDTree(pt, K);
}

template <typename CellType, typename PointType>
CellType &NDTGrid2D<CellType, PointType>::operator[](const PointType &pt)
{
  typename DataGrid::Point point;
  point << pt.x, pt.y;
  return grid_[point];
}

template <typename CellType, typename PointType>
CellType &NDTGrid2D<CellType, PointType>::operator[](const Eigen::Vector2d &pt)
{
  return grid_[pt];
}

template <typename CellType, typename PointType>
const CellType &NDTGrid2D<CellType, PointType>::
operator[](const PointType &pt) const
{
  typename DataGrid::Point point;
  point << pt.x, pt.y;
  return grid_[point];
}

template <typename CellType, typename PointType>
const CellType &NDTGrid2D<CellType, PointType>::
operator[](const Eigen::Vector2d &pt) const
{
  return grid_[pt];
}

template <typename CellType, typename PointType>
bool NDTGrid2D<CellType, PointType>::isInside(const PointType &pt)
{
  typename DataGrid::Point point;
  point << pt.x, pt.y;
  return grid_.cellExists(point);
}

template <typename CellType, typename PointType>
bool NDTGrid2D<CellType, PointType>::isInside(const Eigen::Vector2d &pt) const
{
  return grid_.cellExists(pt);
}

template <typename CellType, typename PointType>
OccupancyGrid NDTGrid2D<CellType, PointType>::createOccupancyGrid() const
{
  OccupancyGrid oc_grid;
  Eigen::Vector2d corner_origin(-static_cast<double>(grid_.left()) * cell_size_,
                                grid_.up() * cell_size_);
  eigt::transform2d_t<double> trans =
      eigt::transBtwFrames(origin_, Eigen::Vector3d(0, 0, 0));
  oc_grid.origin_.head(2) = trans * corner_origin;
  oc_grid.origin_(2) = origin_(2);

  oc_grid.centroid_ << grid_.left() * cell_size_, grid_.up() * cell_size_;

  oc_grid.width_ = grid_.width();
  oc_grid.height_ = grid_.height();
  oc_grid.resolution_ = cell_size_;
  oc_grid.cells_.reserve(oc_grid.width_ * oc_grid.height_);
  for (const auto &cell_it : grid_) {
    if (cell_it.get() == nullptr)
      oc_grid.cells_.push_back(-1);
    else if (cell_it->getOccupancyRaw() > 0.8)
      oc_grid.cells_.push_back(99);
    else
      oc_grid.cells_.push_back(0);
  }
  return oc_grid;
}

template <typename CellType, typename PointType>
typename NDTGrid2D<CellType, PointType>::SelfType
NDTGrid2D<CellType, PointType>::createCoarserGrid(float cell_size) const
{
  if (cell_size < cell_size_)
    cell_size = cell_size_;

  SelfType coarse_grid;
  coarse_grid.setCellSize(cell_size);
  auto cells = grid_.getValidCells();
  coarse_grid.mergeIn(cells, true);
  return coarse_grid;
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
void NDTGrid2D<CellType, PointType>::transform(const Transform &trans)
{
  // get valid cells apply on them transformation and fuse them back to empty
  // grid
  typename DataGrid::CellVector cells = grid_.getValidCells();
  grid_.clear();
  Transform only_rotate = trans;
  only_rotate.matrix().block(0, 2, 2, 1) << 0, 0;
  transformNDTCells(cells, only_rotate);
  for (auto &&cell : cells) {
    if (cell.hasGaussian())
      grid_.addCell(cell.getMean().head(2), std::move(cell), true);
  }
  // transform origin of ndt grid
  Transform only_translate = trans;
  only_translate.linear().setIdentity();
  origin_ = eigt::transformPose(origin_, only_translate);
  updateMeansCloud();
  updateKDTree();
}

template <typename CellType, typename PointType>
typename NDTGrid2D<CellType, PointType>::Transform
NDTGrid2D<CellType, PointType>::move(const Transform &transform)
{
  if (std::abs(origin_(2)) > 1e-7) {
    std::cerr << "[NDT_GRID2D]: you may not use move of grid when grid is not "
                 "aligned with Kartesian coordinates (your grid is rotated). "
                 "Please use transform instead. origin: "
              << origin_.transpose() << std::endl;
    return transform;
  }
  Transform out_trans = transform;
  bool change = false;
  // std::cout << "cumul 1: \n" << out_trans.matrix() << std::endl;

  Eigen::Vector2i move(0, 0);  // move in discreat steps of grid
  Eigen::Vector2d continuous_move(0, 0);

  if (std::floor(std::abs(out_trans(0, 2)) / cell_size_) > 0) {
    double res = out_trans(0, 2) / cell_size_;
    if (res < 0)
      move(0) = std::ceil(res);
    else
      move(0) = std::floor(res);
    continuous_move(0) = move(0) * cell_size_;
    // used movement in x direction is discarted
    out_trans(0, 2) -= move(0) * cell_size_;
    change = true;
  }
  if (std::floor(std::abs(out_trans(1, 2)) / cell_size_) > 0) {
    double res = out_trans(1, 2) / cell_size_;
    if (res < 0)
      move(1) = std::ceil(res);
    else
      move(1) = std::floor(res);
    continuous_move(1) = move(1) * cell_size_;
    // used movement in y direction is discarted
    out_trans(1, 2) -= move(1) * cell_size_;
    change = true;
  }

  if (change) {
    // move underling voxel grid
    grid_.translate(-move, true);
    // move origin of this grid
    origin_(0) += continuous_move(0);
    origin_(1) += continuous_move(1);
    // applies transformation to all means of ndt cells
    for (auto &&cell : getGaussianCells()) {
      // e.g. if grid origin moves one cell forward all old cell needs to
      // recalculate
      // their position one step back
      // std::cout << cell->getMean().transpose();
      cell->setMean(cell->getMean() + Eigen::Vector3d(-continuous_move(0),
                                                      -continuous_move(1), 0));
      // std::cout << "after: " << cell->getMean().transpose() << std::endl;
    }
    updateMeansCloud();
    updateKDTree();

    // std::cout << change << std::endl;
    // std::cout << "move: " << move.transpose() << "   "
    //           << continuous_move.transpose() << std::endl;
    // std::cout << "cumul 2: \n" << out_trans.matrix() << std::endl;
    // std::cout << "origin: " << origin_.transpose() << std::endl;
  }

  return out_trans;
}

template <typename CellType, typename PointType>
NDTGridMsg NDTGrid2D<CellType, PointType>::serialize() const
{
  NDTGridMsg msg;
  Eigen::Vector2d corner_origin(-grid_.left() * cell_size_,
                                grid_.up() * cell_size_);
  eigt::transform2d_t<double> trans =
      eigt::transBtwFrames(origin_, Eigen::Vector3d(0, 0, 0));
  msg.origin_.head(2) = trans * corner_origin;
  msg.size_ << grid_.width() * cell_size_, grid_.height() * cell_size_, 0;
  msg.cell_sizes_ << cell_size_, cell_size_, 0;
  for (auto &&cell : grid_.getValidCellsPtr()) {
    msg.cells_.push_back(cell->serialize());
  }
  return msg;
}

template <typename Cell, typename Pt>
std::ostream &operator<<(std::ostream &os, const NDTGrid2D<Cell, Pt> &grid)
{
  size_t i = 0;
  for (typename NDTGrid2D<Cell, Pt>::DataGrid::ConstIterator cell =
           grid.grid_.cbegin();
       cell != grid.grid_.cend(); ++cell) {
    if (cell->get() == nullptr)
      os << std::setw(1) << "";
    else if ((*cell)->hasGaussian())
      os << std::setw(1) << "W";
    else
      // os << std::setw(3) << std::floor((*cell)->getOccupancyRaw() * 10);
      os << std::setw(1) << ".";
    ++i;
    if (i == grid.grid_.width()) {
      os << std::endl;
      i = 0;
    }
  }
  return os;
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
  std::vector<CellType> occupied_cells = localGrid.grid_.getValidCells();
  // for (auto &&cell : occupied_cells) {
  //   //std::cout << "valid cell: " << cell.toString() << std::endl;
  // }
  return localGrid;
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::transformNDTCells(
    std::vector<CellType> &grid, const Transform &transform)
{
  typename CellType::Transform trans(eigt::convertFromTransform(transform));
  for (CellType &cell : grid) {
    cell.transform(trans);
  }
}

template <typename CellType, typename PointType>
typename NDTGrid2D<CellType, PointType>::CellPtrVector
NDTGrid2D<CellType, PointType>::getKNearestNeighborsVoxel(
    const Eigen::Vector2d &pt, size_t K) const
{
  Eigen::Vector3d point(pt(0), pt(1), 0);
  CellPtrVector res;
  // slowly increase search radius and find sufficient number of points
  CellPtrVector radius_cells = getNeighbors(pt, 1);
  size_t radius = 2;
  while (radius_cells.size() < K && radius < grid_.width()) {
    radius_cells = getNeighbors(pt, radius);
    radius = radius * 2;
  }
  // not enough elements in grid
  if (radius_cells.size() < K)
    return radius_cells;
  // sort elements based on distance from point from parameter
  std::sort(radius_cells.begin(), radius_cells.end(), [&point](CellType *a,
                                                               CellType *b) {
    return (a->getMean() - point).norm() < (a->getMean() - point).norm();
  });
  res.reserve(K);
  std::move(radius_cells.begin(), radius_cells.begin() + K,
            std::back_inserter(res));
  return res;
}

template <typename CellType, typename PointType>
typename NDTGrid2D<CellType, PointType>::CellPtrVector
NDTGrid2D<CellType, PointType>::getKNearestNeighborsKDTree(
    const Eigen::Vector2d &pt, size_t K) const
{
  CellPtrVector res;
  res.reserve(K);
  PointType point;
  point.x = pt(0);
  point.y = pt(1);
  point.z = 0;
  std::vector<int> res_indexes;
  std::vector<float> distances;
  kdtree_.nearestKSearch(point, K, res_indexes, distances);
  for (int id : res_indexes) {
    PointType near_pt = means_->at(id);
    res.push_back(grid_.getCellPtr(Eigen::Vector2d(near_pt.x, near_pt.y)));
  }
  return res;
}

template <typename CellType, typename PointType>
void NDTGrid2D<CellType, PointType>::updateMeansCloud()
{
  typename PointCloud::Ptr pcl(new PointCloud());
  for (auto &&cell : getGaussianCells()) {
    PointType pt;
    Eigen::Vector3d mean = cell->getMean();
    pt.x = mean(0);
    pt.y = mean(1);
    pt.z = mean(2);
    pcl->push_back(pt);
  }
  means_ = pcl;
}

}  // end of namespace slamuk

#endif
