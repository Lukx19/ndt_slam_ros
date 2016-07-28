
#include <gtest/gtest.h>
#include <fstream>
#include <string>
#include <vector>

#include <laser_geometry/laser_geometry.h>
#include <sensor_msgs/LaserScan.h>
#include <sensor_msgs/PointCloud.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/registration/ndt.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl_conversions/pcl_conversions.h>
#include <boost/thread/thread.hpp>

#include <ndt_gslam/registration/ndt2d.h>

#include <Eigen/Dense>

#include <ndt_gslam/utils/eigen_tools.h>

using namespace pcl;

typedef pcl::PointCloud<pcl::PointXYZ> pcl_t;

double EPSILON = 0.001;
double MIN_DISPLACEMENT = 0.4;
double MIN_ROTATION = 0.3;
std::vector<pcl_t::Ptr> scans;
std::vector<Eigen::Vector3d> real_poses;
// std::vector<MatchResult> matches;
NormalDistributionsTransform2DEx<pcl::PointXYZ, pcl::PointXYZ> *matcher;
NormalDistributionsTransform<pcl::PointXYZ, pcl::PointXYZ> *proofer;

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

void testMatch(size_t source_id, size_t target_id)
{
  eigt::pose2d_t<double> real_trans = eigt::getPoseFromTransform(
      eigt::transBtwPoses(real_poses[target_id], real_poses[source_id]));
  // real_trans<<real_trans(0),real_trans(1)+0.077,real_trans(2);
  matcher->setInputSource(scans[source_id]);
  matcher->setInputTarget(scans[target_id]);
  proofer->setInputSource(scans[source_id]);
  proofer->setInputTarget(scans[target_id]);
  pcl_t::Ptr output_pr(new pcl_t);
  pcl_t::Ptr output_m(new pcl_t);

  // Eigen::Matrix4f guess = Eigen::Matrix4f::Identity();
  Eigen::Matrix4f guess =
      eigt::convertFromTransform(
          eigt::transBtwPoses(real_poses[target_id], real_poses[source_id]))
          .cast<float>();
  proofer->align(*output_pr, guess);
  matcher->align(*output_m, guess);
  std::cout << "sucess" << matcher->hasConverged() << std::endl;
  ASSERT_TRUE(matcher->hasConverged());
  std::cout << "result score: " << matcher->getFitnessScore() << std::endl;
  eigt::pose2d_t<double> calc_trans =
      eigt::getPoseFromTransform(eigt::convertToTransform<double>(
          matcher->getFinalTransformation().cast<double>()));
  eigt::pose2d_t<double> proof_trans =
      eigt::getPoseFromTransform(eigt::convertToTransform<double>(
          proofer->getFinalTransformation().cast<double>()));
  std::cout << "PROOF TRANSFORM:" << proof_trans.transpose() << std::endl;
  std::cout << "CALC TRANSFORM:" << calc_trans.transpose() << std::endl;
  eigt::pose2d_t<double> diff = (proof_trans - calc_trans).cwiseAbs();
  std::cout << "DIFF:" << diff.transpose() << std::endl << std::endl;
  EXPECT_LT(diff(0), 0.05);
  EXPECT_LT(diff(1), 0.05);
  EXPECT_LT(diff(2), 0.2);
}

TEST(NormalDistributionsTransform2DEx, scanPairCalcError)
{
  // matches.clear();
  matcher->setResolution(1);
  size_t target = 0;
  std::cout << "start test" << std::endl;
  for (size_t i = 0; i < 499; ++i) {
    eigt::transform2d_t<double> real_trans =
        eigt::transBtwPoses(real_poses[target], real_poses[i]);
    if (eigt::getAngle(real_trans) > MIN_ROTATION ||
        eigt::getDisplacement(real_trans) > MIN_DISPLACEMENT) {
      std::cout << "matching " << i << "  " << target << std::endl;
      testMatch(i, target);
      target = i;
      // return;
    }
  }
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  preparePoseData("data/");
  prepareLaserScans("data/");
  matcher =
      new NormalDistributionsTransform2DEx<pcl::PointXYZ, pcl::PointXYZ>();
  proofer = new NormalDistributionsTransform<pcl::PointXYZ, pcl::PointXYZ>();
  int res = RUN_ALL_TESTS();
  delete matcher;
  delete proofer;
  return res;
}
