#include <graph_slam_uk/graph_slam2d_node.h>
#include <graph_slam_uk/slam_optimizer/ndt_scanmatcher.h>
//#include <graph_slam_uk/gauss_newton_optimalizer2d.h>
#include <graph_slam_uk/slam_optimizer/slam2d.h>

using namespace slamuk;

int main(int argc, char **argv)
{
  ros::init(argc, argv, "graph_slam");
  ros::NodeHandle n;
  ros::NodeHandle n_private("~");
  NdtScanmatcher scanmatcher;
  Slam2D<pcl::PointCloud<pcl::PointXYZ>::Ptr> opt_engine(scanmatcher);

  GraphSlamNode<pcl::PointCloud<pcl::PointXYZ>::Ptr> slam(
      n, n_private, opt_engine, scanmatcher);
  slam.start();

  return 0;
}
