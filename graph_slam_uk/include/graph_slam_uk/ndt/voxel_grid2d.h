#ifndef GRAPH_SLAM_UK_VOXEL_GRID2D
#define GRAPH_SLAM_UK_VOXEL_GRID2D

#include <vector>
#include <Eigen/Dense>
#include <memory>
#include <iostream>
#include <iomanip>

/*
* CellType should have these methods implemented:
* operator+=
  operator=
*/

namespace slamuk
{
template <typename CellType>
class VoxelGrid2D
{
  typedef std::vector<std::unique_ptr<CellType>> CellDataVector;

public:
  typedef typename CellDataVector::iterator Iterator;
  typedef typename CellDataVector::const_iterator ConstIterator;
  typedef std::vector<CellType*> CellPtrVector;
  typedef std::vector<CellType> CellVector;
  typedef Eigen::Vector2d Point;
  typedef Eigen::Matrix<size_t, 2, 1> CentroidIds;
  // typedef std::unique_ptr<CellType> CellTypePtr;
public:
  explicit VoxelGrid2D();

  VoxelGrid2D(const VoxelGrid2D& other);
  VoxelGrid2D(VoxelGrid2D&& other) = default;
  VoxelGrid2D& operator=(const VoxelGrid2D& other);
  VoxelGrid2D& operator=(VoxelGrid2D&& other) = default;
  VoxelGrid2D clone() const;
  // bool addCell(size_t row, size_t col,const CellTypePtr & cell, bool
  // destructive = true);
  // bool addCell(size_t row, size_t col,CellTypePtr && cell, bool destructive =
  // true);
  void addCell(const Point& posistion, const CellType& cell,
               bool merging = true);
  void addCell(const Point& posistion, CellType&& cell, bool merging = true);
  bool removeCell(const Point& posistion);
  // CellType & getCell(size_t row, size_t col);
  bool isInside(const Point& pt) const;

  // CellVector getNeighbors(size_t row, size_t col, size_t radius) const;
  // return all initialized neighbors in radius including cell where pt belongs
  CellPtrVector getNeighbors(const Point& pt, size_t radius) const;

  bool cellExists(const Point& pt) const;
  CellType& operator[](const Point& pt);
  CellType* getCellPtr(const Point& pt) const;
  // returns all initialized cells by values. Should be used for creation of new
  // grids
  CellVector getValidCells() const;
  CellPtrVector getValidCellsPtr() const;
  // enlarges grid in all directions to size based on parameters. If parameters
  // are smaller than current state no enlarging or resizing is done.
  void enlarge(float minx, float miny, float maxx, float maxy)
  {
    enlargeGrid(calcIncLeft(minx), calcIncRight(maxx), calcIncUp(maxy),
                calcIncDown(miny));
  }
  // return all cells laying on the ray. Uninitialized cells are created with
  // default constructor
  CellPtrVector rayTrace(const Point& start, const Point& end);
  // starts at 0,0
  CellPtrVector rayTrace(const Point& end)
  {
    return rayTrace(Point(0, 0), end);
  }
  /**
   * @brief      Translates grid
   *
   * @param[in]  translation  Translation of grid [x,y]. [1,1] translate one
   *cell
   *                          right and one cell up [-1,-1] translate one cell
   *                          left.
   * @param[in]  discard      True will not resize grid
   */
  void translate(const Eigen::Vector2i& translation, bool discard = true);

  /**
   * @brief      Removes all valid cells from structure. All pointers, iterators
   *             are invalidated. Allocated capacity of the grid stays the same.
   */
  void clear();

  /////////////////Parameters
  void setCellSize(float cell_size)
  {
    cell_size_ = cell_size;
    cell_size_half_ = cell_size / 2;
  }

  float getCellSize() const
  {
    return cell_size_;
  }

  size_t width() const
  {
    return width_left_ + width_right_ + 1;
  }
  size_t height() const
  {
    return height_up_ + height_down_ + 1;
  }

  // //////////////////iterators
  Iterator begin()
  {
    return cells_.begin();
  }

  Iterator end()
  {
    return cells_.end();
  }

  ConstIterator cbegin() const
  {
    return cells_.cbegin();
  }

  ConstIterator cend() const
  {
    return cells_.cend();
  }
  size_t validCells() const
  {
    return valid_cells_;
  }
  // output
  template <typename Cell>
  friend std::ostream& operator<<(std::ostream& os,
                                  const VoxelGrid2D<Cell>& grid);

private:
  float cell_size_;
  float cell_size_half_;
  // number of cells in all directions from centroid cell
  size_t width_left_, width_right_, height_up_, height_down_;
  size_t valid_cells_;
  CellDataVector cells_;

  size_t calcIndex(const Point& pt) const;
  size_t calcIndex(const std::pair<size_t, size_t>& coords) const;
  std::pair<size_t, size_t> calcCoordinates(const Point& pt) const;

  void enlargeToFit(const Point& pt);
  size_t calcIncLeft(float pt) const;
  size_t calcIncRight(float pt) const;
  size_t calcIncUp(float pt) const;
  size_t calcIncDown(float pt) const;
  void enlargeGrid(size_t left, size_t right, size_t up, size_t down);
  void moveLeft(size_t move);
  void moveRight(size_t move);
  void moveUp(size_t move);
  void moveDown(size_t move);
};

////////////////////////IMPLEMENTATION//////////
template <typename CellType>
VoxelGrid2D<CellType>::VoxelGrid2D()
  : cell_size_(0.25f)
  , cell_size_half_(0.125f)
  , width_left_(0)
  , width_right_(0)
  , height_up_(0)
  , height_down_(0)
  , valid_cells_(1)
{
  std::unique_ptr<CellType> cell(new CellType());
  cells_.push_back(std::move(cell));
}

template <typename CellType>
VoxelGrid2D<CellType>::VoxelGrid2D(const VoxelGrid2D& other)
  : cell_size_(other.cell_size_)
  , cell_size_half_(other.cell_size_half_)
  , width_left_(other.width_left_)
  , width_right_(other.width_right_)
  , height_up_(other.height_up_)
  , height_down_(other.height_down_)
{
  cells_.resize(other.cells_.size());
  for (const std::unique_ptr<CellType>& cell : other.cells_) {
    CellType* new_cell = new CellType;
    *new_cell = *cell;
    cells_.emplace_back(new_cell);
  }
}

template <typename CellType>
VoxelGrid2D<CellType>& VoxelGrid2D<CellType>::operator=(const VoxelGrid2D& other)
{
  if (&other == this)
    return *this;
  return operator=(VoxelGrid2D<CellType>(other));  // using move constructor
}

template <typename CellType>
VoxelGrid2D<CellType> VoxelGrid2D<CellType>::clone() const
{
  VoxelGrid2D new_grid;
  new_grid.setCellSize(cell_size_);
  return new_grid;
}

template <typename CellType>
void VoxelGrid2D<CellType>::addCell(const Point& position, const CellType& cell,
                                    bool merging)
{
  // resizing to fit new cell
  enlargeToFit(position);
  size_t idx = calcIndex(position);
  if (cells_[idx].get() == nullptr) {
    CellType* new_cell = new CellType;
    *new_cell = cell;
    cells_[idx].reset(new_cell);
    ++valid_cells_;
  } else {
    if (merging) {
      *cells_[idx] += cell;
    } else {
      *cells_[idx] = cell;
    }
  }
}

template <typename CellType>
void VoxelGrid2D<CellType>::addCell(const Point& position, CellType&& cell,
                                    bool merging)
{
  // resizing to fit new cell
  enlargeToFit(position);

  size_t idx = calcIndex(position);
  if (cells_[idx].get() == nullptr) {
    CellType* new_cell = new CellType;
    *new_cell = std::move(cell);
    cells_[idx].reset(new_cell);
    ++valid_cells_;
  } else {
    if (merging) {
      *cells_[idx] += cell;
    } else {
      *cells_[idx] = std::move(cell);
    }
  }
}

template <typename CellType>
bool VoxelGrid2D<CellType>::removeCell(const Point& posistion)
{
  if (!isInside(posistion))
    return false;
  size_t idx = calcIndex(posistion);
  if (cells_[idx].get() != nullptr) {
    --valid_cells_;
  }
  cells_[idx].reset(nullptr);
}

template <typename CellType>
bool VoxelGrid2D<CellType>::isInside(const Point& pt) const
{
  if (pt(0) >= -width_left_ * cell_size_ - cell_size_half_ &&
      pt(0) <= width_right_ * cell_size_ + cell_size_half_ &&
      pt(1) >= -height_down_ * cell_size_ - cell_size_half_ &&
      pt(1) <= height_up_ * cell_size_ + cell_size_half_)
    return true;
  else
    return false;
}

template <typename CellType>
typename VoxelGrid2D<CellType>::CellPtrVector
VoxelGrid2D<CellType>::getNeighbors(const Point& pt, size_t radius) const
{
  typedef std::pair<size_t, size_t> Coordinates;
  CellPtrVector res;
  if (!isInside(pt))
    return res;
  if (radius == 0) {
    if (cellExists(pt))
      res.push_back(cells_[calcIndex(pt)].get());
    return res;
  }
  Coordinates win_center = calcCoordinates(pt);
  long left = static_cast<long>(win_center.first) - static_cast<long>(radius);
  if (left < 0)
    left = 0;
  long up = static_cast<long>(win_center.second) - static_cast<long>(radius);
  if (up < 0)
    up = 0;
  size_t right = win_center.first + radius + 1;
  if (right > width())
    right = width();
  size_t down = win_center.second + radius + 1;
  if (down > height())
    down = height();
  size_t left_unsig = static_cast<size_t>(left);
  size_t up_unsig = static_cast<size_t>(up);
  for (size_t row = up_unsig; row < down; ++row) {
    for (size_t col = left_unsig; col < right; ++col) {
      size_t idx = calcIndex(std::make_pair(col, row));
      if (cells_[idx].get() != nullptr) {
        res.push_back(cells_[idx].get());
      }
    }
  }
  return res;
}

template <typename CellType>
bool VoxelGrid2D<CellType>::cellExists(const Point& pt) const
{
  if (!isInside(pt))
    return false;
  if (cells_[calcIndex(pt)].get() == nullptr)
    return false;
  return true;
}

template <typename CellType>
CellType& VoxelGrid2D<CellType>::operator[](const Point& pt)
{
  return *cells_[calcIndex(pt)];
}

template <typename CellType>
CellType* VoxelGrid2D<CellType>::getCellPtr(const Point& pt) const
{
  return cells_[calcIndex(pt)].get();
}

// returns all initialized cells by values. Should be used for creation of new
// grids

template <typename CellType>
typename VoxelGrid2D<CellType>::CellVector
VoxelGrid2D<CellType>::getValidCells() const
{
  CellVector good_cells;
  good_cells.reserve(valid_cells_);
  for (auto& cell : cells_) {
    if (cell.get() != nullptr) {
      good_cells.push_back(*cell);
    }
  }
  return good_cells;
}

template <typename CellType>
typename VoxelGrid2D<CellType>::CellPtrVector
VoxelGrid2D<CellType>::getValidCellsPtr() const
{
  CellPtrVector good_cells;
  good_cells.reserve(valid_cells_);
  for (auto& cell : cells_) {
    if (cell.get() != nullptr) {
      good_cells.push_back(cell.get());
    }
  }
  return good_cells;
}

// return all cells laying on the ray. Uninitialized cells are created with
// default constructor

template <typename CellType>
typename VoxelGrid2D<CellType>::CellPtrVector
VoxelGrid2D<CellType>::rayTrace(const Point& start, const Point& end)
{
  CellPtrVector res;
  double minx = -width_left_ * cell_size_ - cell_size_half_;
  double maxx = width_right_ * cell_size_ + cell_size_half_;
  double miny = -height_down_ * cell_size_ - cell_size_half_;
  double maxy = height_up_ * cell_size_ + cell_size_half_;

  // fitting end point of raytracing to grid
  Point end_in;
  Point start_in;
  double diff_x = end(0) - start(0);
  double diff_y = end(1) - start(1);
  double slope = diff_y / diff_x;
  if (std::isinf(slope)) {
    // vertical line special case
    if (!isInside(start)) {
      start_in(0) = start(0);
      if (start(1) < 0)
        start_in(1) = miny;
      else
        start_in(1) = maxy;
    } else {
      start_in = start;
    }
    if (!isInside(end)) {
      end_in(0) = end(0);
      if (end(1) < 0)
        end_in(1) = miny;
      else
        end_in(1) = maxy;
    } else {
      end_in = end;
    }
  } else {
    // for all other lines
    if (!isInside(start)) {
      if (start(0) < 0)
        start_in(0) = minx;
      else
        start_in(0) = maxx;
      start_in(1) = slope * start_in(0) + start(1);
    } else {
      start_in = start;
    }

    if (!isInside(end)) {
      if (end(0) < 0)
        end_in(0) = minx;
      else
        end_in(0) = maxx;
      end_in(1) = slope * end_in(0) + start(1);
    } else {
      end_in = end;
    }
  }
  // intersection check
  if (end_in(0) > maxx || end_in(0) < minx || end_in(1) > maxy ||
      end_in(1) < miny)
    return res;
  if (start_in(0) > maxx || start_in(0) < minx || start_in(1) > maxy ||
      start_in(1) < miny)
    return res;
  if (std::isinf(slope) && ((start_in(1) <= miny && end_in(1) <= miny) ||
                            (start_in(1) >= maxy && end_in(1) >= maxy)))
    return res;

  // calculation of delta movement along line
  auto idx1 = calcCoordinates(start_in);
  auto idx2 = calcCoordinates(end_in);
  Point delta;
  double length_x = std::abs(static_cast<double>(idx2.first) -
                             static_cast<double>(idx1.first));
  double length_y = std::abs(static_cast<double>(idx2.second) -
                             static_cast<double>(idx1.second));
  delta(0) = (end_in(0) - start_in(0)) / std::max(length_x, length_y);
  delta(1) = (end_in(1) - start_in(1)) / std::max(length_y, length_x);
  delta *= 0.5;
  Point p = start;
  p += delta;
  size_t idx_old = calcIndex(start_in);
  size_t idx_current = calcIndex(p);
  size_t idx_finish = calcIndex(end_in);

  // appending start element
  if (cells_[idx_old].get() == nullptr)
    cells_[idx_old].reset(new CellType());
  res.push_back(cells_[idx_old].get());

  while (idx_current != idx_finish) {
    // if program  is still in the same cell like last time just step forward
    if (idx_current != idx_old) {
      if (cells_[idx_current].get() == nullptr) {
        cells_[idx_current].reset(new CellType());
      }
      res.push_back(cells_[idx_current].get());
    }
    p += delta;
    idx_old = idx_current;
    idx_current = calcIndex(p);
  }
  // appending end point element
  if (cells_[idx_finish].get() == nullptr)
    cells_[idx_finish].reset(new CellType());
  res.push_back(cells_[idx_finish].get());

  return res;
}

// [1,1] translate one cell right and one cell up
// [-1,-1] translate one cell left and one cell down
// discard = false will resize grid

template <typename CellType>
void VoxelGrid2D<CellType>::translate(const Eigen::Vector2i& translation,
                                      bool discard)
{
  size_t left, right, up, down;
  left = right = up = down = 0;
  if (translation(0) < 0)
    left = static_cast<size_t>(std::abs(translation(0)));
  else
    right = static_cast<size_t>(translation(0));

  if (translation(1) < 0)
    down = static_cast<size_t>(translation(1) * -1);
  else
    up = static_cast<size_t>(translation(1));
  if (!discard) {
    enlargeGrid(left, right, up, down);
  }
  moveLeft(left);
  moveRight(right);
  moveUp(up);
  moveDown(down);
}

template <typename CellType>
void VoxelGrid2D<CellType>::clear()
{
  for (auto&& cell : cells_) {
    cell.reset(nullptr);
  }
}

template <typename CellType>
std::ostream& operator<<(std::ostream& os, const VoxelGrid2D<CellType>& grid)
{
  size_t width = grid.width();
  size_t i = 0;
  for (auto&& cell : grid.cells_) {
    if (i == width) {
      os << std::endl;
      i = 0;
    }
    ++i;
    if (cell.get() == nullptr)
      os << std::setw(5) << "NULL";
    else
      os << std::setw(5) << *cell;
  }
  return os;
}

// /////////////////PROTECTED

template <typename CellType>
size_t VoxelGrid2D<CellType>::calcIndex(const Point& pt) const
{
  return calcIndex(calcCoordinates(pt));
}

template <typename CellType>
size_t
VoxelGrid2D<CellType>::calcIndex(const std::pair<size_t, size_t>& coords) const
{
  // row * row_width + coll
  return coords.second * width() + coords.first;
}

template <typename CellType>
std::pair<size_t, size_t>
VoxelGrid2D<CellType>::calcCoordinates(const Point& pt) const
{
  std::pair<size_t, size_t> coords;
  float minx = cell_size_ * width_left_;
  float maxy = cell_size_ * height_up_;
  coords.first = static_cast<size_t>(
      std::floor((pt(0) + minx + cell_size_half_) / cell_size_));
  coords.second = static_cast<size_t>(
      std::floor((-pt(1) + maxy + cell_size_half_) / cell_size_));
  return coords;
}

template <typename CellType>
void VoxelGrid2D<CellType>::enlargeToFit(const Point& pt)
{
  // keep in sync with calculation in isInside. Mostly operator< and widths
  // calculations
  if (isInside(pt))
    return;
  enlargeGrid(calcIncLeft(pt(0)), calcIncRight(pt(0)), calcIncUp(pt(1)),
              calcIncDown(pt(1)));
}

template <typename CellType>
size_t VoxelGrid2D<CellType>::calcIncLeft(float pt) const
{
  float minx = -width_left_ * cell_size_ - cell_size_half_;
  if (pt < minx)
    return static_cast<size_t>(std::ceil((minx - pt) / cell_size_));
  return 0;
}

template <typename CellType>
size_t VoxelGrid2D<CellType>::calcIncRight(float pt) const
{
  float maxx = width_right_ * cell_size_ + cell_size_half_;
  if (pt > maxx)
    return static_cast<size_t>(std::ceil((pt - maxx) / cell_size_));
  return 0;
}

template <typename CellType>
size_t VoxelGrid2D<CellType>::calcIncUp(float pt) const
{
  float maxy = height_up_ * cell_size_ + cell_size_half_;
  if (pt > maxy)
    return static_cast<size_t>(std::ceil((pt - maxy) / cell_size_));
  return 0;
}

template <typename CellType>
size_t VoxelGrid2D<CellType>::calcIncDown(float pt) const
{
  float miny = -height_down_ * cell_size_ - cell_size_half_;
  if (pt < miny)
    return static_cast<size_t>(std::ceil((miny - pt) / cell_size_));
  return 0;
}

// increase in each dimension
template <typename CellType>
void VoxelGrid2D<CellType>::enlargeGrid(size_t left, size_t right, size_t up,
                                        size_t down)
{
  // std::cout << "enlarge: l: " << left << " r: " << right << " u: " << up
  // << " d: " << down << std::endl;
  if (left + right + up + down == 0)
    return;
  CellDataVector large_grid;
  size_t new_size = (width() + left + right) * (height() + up + down);
  large_grid.resize(new_size);
  size_t w = width();
  size_t w_used = 0;
  size_t new_w = width() + left + right;
  size_t offset_x = 0;
  for (size_t i = 0; i < cells_.size(); ++i) {
    size_t new_index = 0;
    if (up > 0 && left == 0) {
      new_index = i + up * (new_w)+offset_x;
    } else if (up == 0 && left == 0) {
      new_index = i + offset_x;
    } else {
      new_index = i + up * (new_w)+left + offset_x;
    }
    large_grid[new_index] = std::move(cells_[i]);
    ++w_used;
    if (w_used == w) {
      offset_x += right + left;
      w_used = 0;
    }
  }
  cells_ = std::move(large_grid);
  height_up_ += up;
  height_down_ += down;
  width_left_ += left;
  width_right_ += right;
}

template <typename CellType>
void VoxelGrid2D<CellType>::moveLeft(size_t move)
{
  if (move == 0)
    return;
  size_t width = this->width();
  if (move > width)
    move = width;
  size_t used, unused;
  used = unused = 0;
  size_t good_cells = width - move;
  for (size_t i = 0; i < cells_.size(); ++i) {
    long idx = static_cast<long>(i) + static_cast<long>(move);
    if (idx < cells_.size() && used < good_cells && unused == 0) {
      ++used;
      unused = 0;
      cells_[i] = std::move(cells_[idx]);
    } else {
      used = 0;
      ++unused;
      if (unused == move)
        unused = 0;
      cells_[i].reset(nullptr);
    }
  }
}

template <typename CellType>
void VoxelGrid2D<CellType>::moveRight(size_t move)
{
  if (move == 0)
    return;
  size_t width = this->width();
  if (move > width)
    move = width;
  size_t used, unused;
  used = unused = 0;
  size_t good_cells = width - move;
  for (long i = cells_.size() - 1; i >= 0; --i) {
    long idx = i - static_cast<long>(move);
    if (idx >= 0 && used < good_cells && unused == 0) {
      ++used;
      unused = 0;
      cells_[i] = std::move(cells_[idx]);
    } else {
      used = 0;
      ++unused;
      if (unused == move)
        unused = 0;
      cells_[i].reset(nullptr);
    }
  }
}

template <typename CellType>
void VoxelGrid2D<CellType>::moveUp(size_t move)
{
  if (move == 0)
    return;
  size_t width = this->width();
  if (move > this->height()) {
    move = this->height();
  }
  for (size_t i = 0; i < cells_.size(); ++i) {
    size_t idx = i + move * width;
    if (idx < cells_.size()) {
      cells_[i] = std::move(cells_[idx]);
    } else {
      cells_[i].reset(nullptr);
    }
  }
}

template <typename CellType>
void VoxelGrid2D<CellType>::moveDown(size_t move)
{
  if (move == 0)
    return;
  size_t width = this->width();
  if (move > this->height()) {
    move = this->height();
  }
  for (long i = cells_.size() - 1; i >= 0; --i) {
    long idx = static_cast<long>(i) -
               static_cast<long>(move) * static_cast<long>(width);
    if (idx >= 0) {
      cells_[i] = std::move(cells_[idx]);
    } else {
      cells_[i].reset(nullptr);
    }
  }
}

}  // end of slamuk namespace
#endif
