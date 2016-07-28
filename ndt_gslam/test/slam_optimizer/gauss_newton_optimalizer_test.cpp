#include <fstream>
#include <string>
#include <vector>

#include <laser_geometry/laser_geometry.h>
#include <sensor_msgs/LaserScan.h>
#include <sensor_msgs/PointCloud.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>

#include <Eigen/Dense>

#include <ndt_gslam/slam_optimizer/gauss_newton_optimalizer2d.h>
#include <ndt_gslam/slam_optimizer/ndt_scanmatcher.h>
#include <ndt_gslam/utils/eigen_tools.h>

using namespace pcl;

typedef pcl::PointCloud<pcl::PointXYZ> pcl_t;

std::vector<pcl_t::Ptr> scans;
std::vector<Eigen::Vector3d> real_poses;
size_t idx = 0;
slamuk::NdtScanmatcher matcher;
slamuk::GaussNewtonOptimalize2d<pcl_t::Ptr> graph(matcher);

std::vector<std::string> split(std::string data, std::string token)
{
  std::vector<std::string> output;
  size_t pos = std::string::npos;
  do {
    pos = data.find(token);
    output.push_back(data.substr(0, pos));
    if (std::string::npos != pos)
      data = data.substr(pos + token.size());
  } while (std::string::npos != pos);
  return output;
}

void prepareLaserScans(std::string folder)
{
  sensor_msgs::PointCloud2 laser_pcl_msg;
  pcl_t laser_pcl;
  laser_geometry::LaserProjection projector;

  scans.clear();
  std::fstream laser_msg;
  std::stringstream file;
  file << folder << "/data.scans";
  laser_msg.open(file.str().c_str(), std::ios_base::in);
  if (!laser_msg) {
    std::cout << "File" << folder << "/data.scans NOT FOUND" << std::endl;
  }
  std::string line;
  unsigned int seq = 0;
  while (std::getline(laser_msg, line)) {
    std::vector<std::string> parts = split(line, ",");
    sensor_msgs::LaserScan msg;
    msg.header.seq = seq;
    msg.header.frame_id = "base_link";
    msg.angle_min = static_cast<float>(std::atof(parts[1].c_str()));
    msg.angle_max = 2.26456475258f;
    msg.range_min = 0.0230000000447f;
    msg.range_max = 60.0f;
    msg.time_increment = 1.73611115315e-05f;
    msg.angle_increment = static_cast<float>(std::atof(parts[2].c_str()));
    for (size_t i = 4; i < parts.size(); ++i) {
      msg.ranges.push_back(static_cast<float>(std::atof(parts[i].c_str())));
    }
    // msg to pcl
    projector.projectLaser(msg, laser_pcl_msg);
    pcl::moveFromROSMsg(laser_pcl_msg, laser_pcl);
    scans.push_back(laser_pcl.makeShared());
    ++seq;
  }
  laser_msg.close();
}

void preparePoseData(std::string folder)
{
  real_poses.clear();
  std::fstream pose_msg;
  std::stringstream file;
  file << folder << "/data.poses";
  pose_msg.open(file.str().c_str(), std::ios_base::in);
  if (!pose_msg) {
    std::cout << "File" << folder << "/data.poses NOT FOUND" << std::endl;
  }
  std::string line;
  Eigen::Vector3d pose;
  while (std::getline(pose_msg, line)) {
    std::vector<std::string> parts = split(line, ",");
    pose << std::atof(parts[1].c_str()), std::atof(parts[2].c_str()),
        std::atof(parts[3].c_str());
    real_poses.push_back(pose);
  }
  pose_msg.close();
}

void test()
{
  graph.addPose(real_poses[0], scans[0]);
  for (size_t i = 1; i < real_poses.size(); ++i) {
    auto real_trans = eigt::transBtwPoses(real_poses[i - 1], real_poses[i]);
    slamuk::MatchResult res;
    // res =matcher.match(scans[i], scans[i-1],real_trans.matrix());
    graph.addPose(real_poses[i], scans[i]);
    // std::cout<<res.inform_.inverse()<<std::endl;
    graph.addLastConstrain(eigt::getPoseFromTransform(real_trans),
                           Eigen::Matrix3d::Identity());
    std::cout << "Added nodes: " << i - 1 << ":" << i << std::endl;
    std::cout << "Real displacement: " << eigt::getDisplacement(real_trans)
              << std::endl;
    graph.tryLoopClose();
    graph.getGraphSerialized(std::cout);
  }
}

int main(int argc, char **argv)
{
  std::vector<std::string> args(argv, argv + argc);
  if (args.size() != 2) {
    std::cout << "Correct format of arguments: \n Path to the folder with "
                 "data.poses, data.scans was not provided";
    std::cout << std::endl;
    return 0;
  }
  preparePoseData(args[1]);
  prepareLaserScans(args[1]);
  test();
  return 0;
}