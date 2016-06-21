#ifndef POINT_CLOUD_TOOLS
#define POINT_CLOUD_TOOLS

#include <pcl/point_cloud.h>

namespace pcl
{
template <typename PointType, typename Scalar>
void getMinMax2D(const pcl::PointCloud<PointType> &pcl, Scalar &minx,
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
}

#endif