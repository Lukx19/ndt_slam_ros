#include <gtest/gtest.h>
#include <ndt_gslam/ndt/cell_policy2d.h>
#include <ndt_gslam/ndt/ndt_cell.h>
#include <ndt_gslam/ndt/voxel_grid2d.h>

#include <iostream>

using namespace slamuk;

typedef VoxelGrid2D<int> GridSimple;
typedef VoxelGrid2D<NDTCell<CellPolicy2d>> GridCell;
typedef Eigen::Vector2d Point;

GridSimple createGridSimple()
{
  GridSimple grid;
  grid.setCellSize(1);
  int i = 0;
  for (double pty = 10; pty > -4; pty -= 1) {
    for (double ptx = -5; ptx < 8; ptx += 1) {
      grid.addCell(Point(ptx, pty), i);
      ++i;
    }
  }
  return grid;
}
TEST(VoxelGrid2D, enlarge)
{
  GridSimple grid;
  grid.setCellSize(1);
  grid.enlarge(0, 0, 10, 10);
  EXPECT_EQ(11, grid.width());
  EXPECT_EQ(11, grid.height());
  grid.enlarge(-10, -10, 5, 5);
  EXPECT_EQ(21, grid.width());
  EXPECT_EQ(21, grid.height());
  grid.enlarge(-10, -10, 11, 5);
  EXPECT_EQ(22, grid.width());
  EXPECT_EQ(21, grid.height());
  EXPECT_TRUE(grid.isInside(Point(10, 10)));
  grid.enlarge(100, 100, -100, -100);
  EXPECT_EQ(22, grid.width());
  EXPECT_EQ(21, grid.height());

  grid = createGridSimple();
  int zero = grid[Point(0, 0)];
  // -6 -5 8 11
  grid.enlarge(-6, -4, 8, 11);
  EXPECT_EQ(zero, grid[Point(0, 0)]);

  // std::cout << grid << std::endl << std::endl;

  // grid = createGridSimple();
  // // -6 -5 8 11
  // grid.enlarge(-6, 0, 0, 0);
  // std::cout << grid << std::endl << std::endl;

  // grid = createGridSimple();
  // // -6 -5 8 11
  // grid.enlarge(0, -4, 0, 0);
  // std::cout << grid << std::endl << std::endl;

  // grid = createGridSimple();
  // // -6 -5 8 11
  // grid.enlarge(0, 0, 8, 0);
  // std::cout << grid << std::endl << std::endl;

  // grid = createGridSimple();
  // // -6 -5 8 11
  // grid.enlarge(0, 0, 0, 11);
  // std::cout << grid << std::endl << std::endl;

  // grid = createGridSimple();
  // // -6 -5 8 11
  // grid.enlarge(-6, -4, 0, 0);
  // std::cout << grid << std::endl << std::endl;

  // grid = createGridSimple();
  // // -6 -5 8 11
  // grid.enlarge(-6, 0, 8, 0);
  // std::cout << grid << std::endl << std::endl;

  // grid = createGridSimple();
  // // -6 -5 8 11
  // grid.enlarge(-6, 0, 0, 11);
  // std::cout << grid << std::endl << std::endl;

  // grid = createGridSimple();
  // // -6 -5 8 11
  // grid.enlarge(0, -4, 0, 11);
  // std::cout << grid << std::endl << std::endl;

  // grid = createGridSimple();
  // // -6 -5 8 11
  // grid.enlarge(0, -4, 8, 0);
  // std::cout << grid << std::endl << std::endl;

  // grid = createGridSimple();
  // // -6 -5 8 11
  // grid.enlarge(0, 0, 8, 11);
  // std::cout << grid << std::endl << std::endl;
}

TEST(VoxelGrid2D, cellGetSetValidate)
{
  GridSimple grid;
  grid.setCellSize(1);
  grid.addCell(Point(0, 0), 10, false);
  EXPECT_EQ(10, grid[Point(0, 0)]);
  grid.addCell(Point(0, 0), 2, true);
  EXPECT_EQ(12, grid[Point(0, 0)]);
  // lavalue interface
  grid.addCell(Point(0, 0), 0, false);
  int val = 10;
  int val2 = 2;
  grid.addCell(Point(0, 0), val, false);
  EXPECT_EQ(10, grid[Point(0, 0)]);
  grid.addCell(Point(0, 0), val2, true);
  EXPECT_EQ(12, grid[Point(0, 0)]);
  std::cout << grid << std::endl;
  // adding cell to other parts of grid than origin
  grid.addCell(Point(2, 2), 5, false);
  std::cout << grid << std::endl;
  EXPECT_EQ(3, grid.width());
  EXPECT_EQ(3, grid.height());
  EXPECT_EQ(5, grid[Point(2, 2)]);
  grid.addCell(Point(-2, -2), -5, false);
  std::cout << grid << std::endl;
  EXPECT_EQ(5, grid.width());
  EXPECT_EQ(-5, grid[Point(-2, -2)]);
  EXPECT_EQ(5, grid[Point(2, 2)]);
  grid.addCell(Point(-4, 2), 5, false);
  grid.addCell(Point(-4, 4), 4, false);
  EXPECT_EQ(7, grid.width());
  EXPECT_EQ(7, grid.height());
  EXPECT_TRUE(grid.cellExists(Point(-3.8, 2.2)));
  EXPECT_TRUE(grid.isInside(Point(-4.4999, 2.49999)));
  EXPECT_FALSE(grid.cellExists(Point(-3.4, 2.6)));
  EXPECT_TRUE(grid.isInside(Point(-1.6999, 2.49999)));
  EXPECT_FALSE(grid.isInside(Point(-8.4999, 2.49999)));
  EXPECT_FALSE(grid.cellExists(Point(-8.4999, 2.49999)));
  std::cout << grid << std::endl;
  std::vector<int> valid_cells = {4, 5, 5, 12, -5};
  size_t i = 0;
  for (auto &&cell : grid.getValidCells()) {
    EXPECT_EQ(valid_cells[i], cell);
    ++i;
  }
  i = 0;
  for (auto &&cell : grid.getValidCellsPtr()) {
    EXPECT_EQ(valid_cells[i], *cell);
    ++i;
  }
  EXPECT_EQ(valid_cells.size(), grid.validCells());
}

TEST(VoxelGrid2D, translation)
{
  GridSimple grid = createGridSimple();
  int zero = grid[Point(0, 0)];
  std::cout << grid << std::endl << std::endl;
  // translate up and discard
  grid.translate(Eigen::Vector2i(0, 1), true);
  EXPECT_EQ(13, grid[Point(-5, 10)]);
  grid = createGridSimple();
  grid.translate(Eigen::Vector2i(0, 13), true);
  EXPECT_EQ(169, grid[Point(-5, 10)]);
  grid = createGridSimple();
  grid.translate(Eigen::Vector2i(0, 20), true);
  EXPECT_FALSE(grid.cellExists(Point(-5, 10)));

  // translate down and discard
  grid = createGridSimple();
  grid.translate(Eigen::Vector2i(0, -1), true);
  EXPECT_EQ(156, grid[Point(-5, -3)]);
  grid = createGridSimple();
  grid.translate(Eigen::Vector2i(0, -13), true);
  EXPECT_EQ(0, grid[Point(-5, -3)]);
  grid = createGridSimple();
  grid.translate(Eigen::Vector2i(0, -20), true);
  EXPECT_FALSE(grid.cellExists(Point(-5, -3)));
  grid.translate(Eigen::Vector2i(0, -1), true);
  grid = createGridSimple();

  // // non discarding
  // grid = createGridSimple();
  // grid.translate(Eigen::Vector2i(1, -3), true);
  // std::cout << grid << std::endl << std::endl;

  // grid = createGridSimple();
  // grid.translate(Eigen::Vector2i(1, -3), false);
  // std::cout << grid << std::endl << std::endl;
}

TEST(VoxelGrid2D, getNeighbours)
{
  GridSimple grid;
  grid.setCellSize(1);
  grid.addCell(Point(0, 0), 0);
  grid.addCell(Point(-1, -1), -111);
  grid.addCell(Point(1, -1), 111);
  grid.addCell(Point(-1, 1), -101);
  grid.addCell(Point(1, 1), 101);

  grid.addCell(Point(-2, -2), -222);
  grid.addCell(Point(2, -2), 222);
  grid.addCell(Point(-2, 2), -202);
  grid.addCell(Point(2, 2), 202);

  grid.addCell(Point(-3, -3), -333);
  grid.addCell(Point(3, -3), 333);
  grid.addCell(Point(-3, 3), -303);
  grid.addCell(Point(3, 3), 303);
  std::cout << grid << std::endl << std::endl;
  std::vector<int> valid_res = {0};
  size_t i = 0;
  for (auto &&ng : grid.getNeighbors(Point(0, 0), 0)) {
    EXPECT_EQ(valid_res[i], *ng);
    ++i;
  }
  i = 0;
  valid_res = {-101, 101, 0, -111, 111};
  for (auto &&ng : grid.getNeighbors(Point(0, 0), 1)) {
    EXPECT_EQ(valid_res[i], *ng);
    ++i;
  }
  i = 0;
  valid_res = {-202, 202, -101, 101, 0, -111, 111, -222, 222};
  for (auto &&ng : grid.getNeighbors(Point(0, 0), 2)) {
    EXPECT_EQ(valid_res[i], *ng);
    ++i;
  }
  i = 0;
  valid_res = {-303, 303, -202, 202, -101, 101, 0,
               -111, 111, -222, 222, -333, 333};
  for (auto &&ng : grid.getNeighbors(Point(0, 0), 10)) {
    EXPECT_EQ(valid_res[i], *ng);
    ++i;
  }
  i = 0;
  valid_res = {303, 202};
  for (auto &&ng : grid.getNeighbors(Point(3, 3), 1)) {
    EXPECT_EQ(valid_res[i], *ng);
    ++i;
  }
  i = 0;
  valid_res = {303, 202, 101};
  for (auto &&ng : grid.getNeighbors(Point(3, 3), 2)) {
    EXPECT_EQ(valid_res[i], *ng);
    ++i;
  }
  i = 0;
  valid_res = {303, 202, 101, 0};
  for (auto &&ng : grid.getNeighbors(Point(2, 2), 2)) {
    EXPECT_EQ(valid_res[i], *ng);
    ++i;
  }
  i = 0;

  valid_res = {0, -111, -222, -333};
  for (auto &&ng : grid.getNeighbors(Point(-2, -2), 2)) {
    EXPECT_EQ(valid_res[i], *ng);
    ++i;
  }
  i = 0;
}

TEST(VoxelGrid2D, rayTracing)
{
  GridSimple grid = createGridSimple();
  grid.setCellSize(1);
  grid.addCell(Point(0, 0), 0, false);
  std::cout << grid << std::endl << std::endl;
  std::vector<int> valid_res = {0, 136, 123, 124, 111, 112, 99};
  size_t i = 0;
  for (auto &&cell : grid.rayTrace(Point(3, 3))) {
    EXPECT_EQ(valid_res[i], *cell);
    ++i;
  }

  valid_res = {0, 134, 121, 120, 119};
  i = 0;
  for (auto &&cell : grid.rayTrace(Point(-3.2, 1.4))) {
    EXPECT_EQ(valid_res[i], *cell);
    ++i;
  }

  valid_res = {0, 134, 133, 132, 131, 130};
  i = 0;
  for (auto &&cell : grid.rayTrace(Point(-90, -1.3))) {
    EXPECT_EQ(valid_res[i], *cell);
    ++i;
  }

  valid_res = {0, 122, 109, 96, 83, 70, 57, 44, 31, 18, 5};
  i = 0;
  for (auto &&cell : grid.rayTrace(Point(0, 20))) {
    EXPECT_EQ(valid_res[i], *cell);
    ++i;
  }

  EXPECT_EQ(0, grid.rayTrace(Point(3, -50), Point(3, -30)).size());
  EXPECT_EQ(0, grid.rayTrace(Point(3, -30), Point(3, -50)).size());
  EXPECT_EQ(0, grid.rayTrace(Point(3, 20), Point(3, 30)).size());
  EXPECT_EQ(0, grid.rayTrace(Point(3, 30), Point(3, 20)).size());
  EXPECT_EQ(0, grid.rayTrace(Point(-200, 0), Point(-100, 0)).size());
  EXPECT_EQ(0, grid.rayTrace(Point(-100, 0), Point(-200, 0)).size());
  EXPECT_EQ(0, grid.rayTrace(Point(200, 0), Point(100, 0)).size());
  EXPECT_EQ(0, grid.rayTrace(Point(100, 0), Point(200, 0)).size());

  // diagonal out of scope rays
  EXPECT_EQ(0, grid.rayTrace(Point(20, 0), Point(30, 30)).size());
  EXPECT_EQ(0, grid.rayTrace(Point(30, 30), Point(20, 0)).size());
  EXPECT_EQ(0, grid.rayTrace(Point(20, 30), Point(30, 0)).size());
  EXPECT_EQ(0, grid.rayTrace(Point(30, 0), Point(20, 30)).size());

  GridSimple grid2;
  grid2.setCellSize(1);
  grid2.addCell(Point(0, 0), 0, false);
  grid2.addCell(Point(-1, -1), -111, false);
  grid2.addCell(Point(1, -1), 111, false);
  grid2.addCell(Point(-1, 1), -101, false);
  grid2.addCell(Point(1, 1), 101, false);

  grid2.addCell(Point(-2, -2), -222, false);
  grid2.addCell(Point(2, -2), 222, false);
  grid2.addCell(Point(-2, 2), -202, false);
  grid2.addCell(Point(2, 2), 202, false);

  grid2.addCell(Point(-3, -3), -333, false);
  grid2.addCell(Point(3, -3), 333, false);
  grid2.addCell(Point(-3, 3), -303, false);
  grid2.addCell(Point(3, 3), 303, false);
  std::cout << grid2 << std::endl << std::endl;
  EXPECT_EQ(4, grid2.rayTrace(Point(-3, 0)).size());
  EXPECT_EQ(4, grid2.rayTrace(Point(-3, 0), Point(-3, 10)).size());
  EXPECT_EQ(9, grid2.rayTrace(Point(-3, -2), Point(3, 2)).size());
  EXPECT_EQ(5, grid2.rayTrace(Point(-2, -3), Point(2, -3)).size());
  EXPECT_EQ(0, grid2.rayTrace(Point(-8, -3), Point(10, -3)).size());
  EXPECT_EQ(7, grid2.rayTrace(Point(-3, -3), Point(10, -3)).size());
  EXPECT_EQ(4, grid2.rayTrace(Point(0, 10)).size());
  std::cout << grid2 << std::endl << std::endl;
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
