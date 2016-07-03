#ifndef GRAPH_SLAM_UK_POINT_CLOUD_TOOLS
#define GRAPH_SLAM_UK_POINT_CLOUD_TOOLS

#include <pcl/point_cloud.h>
#include <vector>

namespace pcl
{
template <typename PointType, typename Scalar>
void getMinMax2D(const pcl::PointCloud<PointType> &pcl, Scalar &minx,
                 Scalar &miny, Scalar &maxx, Scalar &maxy);

template <typename CellType, typename Scalar>
void getMinMaxNDT2D(const std::vector<CellType> &cells, Scalar &minx,
                    Scalar &miny, Scalar &maxx, Scalar &maxy);

////////////////////IMPLEMENTATION ///////////////
template <typename PointType, typename Scalar>
void getMinMax2D(const pcl::PointCloud<PointType> &pcl, Scalar &minx,
                 Scalar &miny, Scalar &maxx, Scalar &maxy)
{
  minx = miny = maxx = maxy = 0;
  for (size_t i = 0; i < pcl.size(); ++i) {
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

template <typename CellType, typename Scalar>
void getMinMaxNDT2D(const std::vector<CellType> &cells, Scalar &minx,
                    Scalar &miny, Scalar &maxx, Scalar &maxy)
{
  minx = miny = maxx = maxy = 0;
  for (auto &&cell : cells) {
    auto mean = cell.getMean();
    if (mean(0) < minx)
      minx = mean(0);
    if (mean(0) > maxx)
      maxx = mean(0);
    if (mean(1) < miny)
      miny = mean(1);
    if (mean(1) > maxy)
      maxy = mean(1);
  }
}
}

#endif