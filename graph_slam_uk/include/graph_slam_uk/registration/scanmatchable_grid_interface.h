#ifndef GRAPH_SLAM_UK_SCANMATCHABLE_GRID_INTERFACE
#define GRAPH_SLAM_UK_SCANMATCHABLE_GRID_INTERFACE

#include <vector>
#include <Eigen/Dense>

namespace slamuk
{
template <typename CellType>
class IScanmatchableGrid
{
public:
  typedef std::vector<CellType *> CellPtrVector;

  /**
   * @brief      Gets all cells with gaussian dist. which lies in square with
   *             length of side 2 x radius and center in pt.
   *
   * @param[in]  pt      Center of interest.
   * @param[in]  radius  Number of cells from the point outwards horizinatly.
   *
   * @return     The neighbors.
   */
  CellPtrVector getNeighbors(const Eigen::Vector2d &pt,
                             size_t radius) const = 0;
  /**
   * @brief      Creates a coarser grid.
   *
   * @param[in]  multiple  The multiple 2 means 2x2 subgrid in original grid
   *                       will create one field in returned grid. Multiple 4
   *                       means 4x4 etc.
   *
   * @return     Returns new grid with coarser layout.
   */
  IScanmatchableGrid<CellType> createCoarserGrid(size_t multiple) const = 0;

  /**
   * @brief      Gets the centroids of cells with gaussian dist. in it.
   *
   * @return     Vector contains centroids of cells with gaussian distribution
   */
  std::vector<Eigen::Vector2d> getGaussianCentroids() const = 0;
};
}

#endif
