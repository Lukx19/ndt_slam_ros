
#include <chrono>
#include <fstream>
#include <string>
#include <vector>

#include <laser_geometry/laser_geometry.h>
#include <sensor_msgs/LaserScan.h>
#include <sensor_msgs/PointCloud.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/registration/ia_ransac.h>
#include <pcl/registration/icp.h>
#include <pcl/registration/ndt.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl_conversions/pcl_conversions.h>
#include <boost/thread/thread.hpp>

#include <ndt_gslam/registration/correlative_estimation2d.h>
#include <ndt_gslam/registration/d2d_ndt2d.h>
#include <ndt_gslam/registration/d2d_ndt2d_robust.h>
#include <ndt_gslam/registration/ndt2d.h>

#include <ndt_gslam/ndt/cell_policy2d.h>
#include <ndt_gslam/ndt/ndt_cell.h>
#include <ndt_gslam/ndt/ndt_grid2d.h>

#include <Eigen/Dense>

#include <ndt_gslam/utils/eigen_tools.h>
#include <ndt_gslam/utils/point_cloud_tools.h>

using namespace pcl;
using namespace slamuk;

typedef pcl::PointCloud<pcl::PointXYZ> pcl_t;
typedef NDTGrid2D<NDTCell<CellPolicy2d>, pcl::PointXYZ> GridType;
typedef GridType::ConstPtr GridTypeConstPtr;
typedef GridType::Ptr GridTypePtr;
typedef NormalDistributionsTransform2DEx<pcl::PointXYZ, pcl::PointXYZ>
    P2DMatcher;
typedef D2DNormalDistributionsTransform2D<pcl::PointXYZ, pcl::PointXYZ>
    D2DMatcher;
typedef CorrelativeEstimation<pcl::PointXYZ, pcl::PointXYZ> CorrMatcher;
typedef D2DNormalDistributionsTransform2DRobust<pcl::PointXYZ, pcl::PointXYZ>
    RobustMatcher;

double EPSILON = 0.001;
double MIN_DISPLACEMENT = 0.4;
double MIN_ROTATION = 0.3;
size_t MAX_LINES = 1000000;
std::vector<pcl_t::Ptr> scans;
std::vector<Eigen::Vector3d> real_poses;
// std::vector<MatchResult> matches;
Registration<pcl::PointXYZ, pcl::PointXYZ> *matcher;
// D2DNormalDistributionsTransform2D<pcl::PointXYZ, pcl::PointXYZ> *proofer;
IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> *proofer;
// SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ> * proofer;
std::string mode_type;

Eigen::Vector3d total_error_ = Eigen::Vector3d::Identity();

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
  size_t line_count = 0;
  while (std::getline(laser_msg, line) && line_count < MAX_LINES) {
    ++line_count;
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
  size_t line_count = 0;
  while (std::getline(pose_msg, line) && line_count < MAX_LINES) {
    ++line_count;
    std::vector<std::string> parts = split(line, ",");
    pose << std::atof(parts[1].c_str()), std::atof(parts[2].c_str()),
        std::atof(parts[3].c_str());
    real_poses.push_back(pose);
  }
  pose_msg.close();
}

void testMatch(size_t source_id, size_t target_id)
{
  std::chrono::time_point<std::chrono::system_clock> s_proof, e_proof, s_match,
      e_match;
  std::chrono::duration<double> elapsed_seconds;
  pcl_t::Ptr output_pr(new pcl_t);
  pcl_t::Ptr output_m(new pcl_t);

  // prepare initial guess
  Eigen::Matrix4f odometry =
      eigt::convertFromTransform(
          eigt::getTransFromPose(real_poses[target_id]).inverse() *
          eigt::getTransFromPose(real_poses[source_id]))
          .cast<float>();
  Eigen::Vector3d real_trans = eigt::getPoseFromTransform(
      eigt::getTransFromPose(real_poses[target_id]).inverse() *
      eigt::getTransFromPose(real_poses[source_id]));
  ////////////////////////////////////////////////////////
  // calculate my matcher
  s_match = std::chrono::system_clock::now();
  if (mode_type == "proof") {
    matcher->setInputSource(scans[source_id]);
    matcher->setInputTarget(scans[target_id]);
  } else {
    GridTypePtr target_grid(new GridType());
    GridTypePtr source_grid(new GridType());
    target_grid->setCellSize(0.25);
    source_grid->setCellSize(0.25);
    target_grid->initializeSimple(*scans[target_id]);
    source_grid->initializeSimple(*scans[source_id]);
    if (mode_type == "d2d") {
      static_cast<D2DMatcher *>(matcher)->setOulierRatio(0.4);
      static_cast<D2DMatcher *>(matcher)->setStepSize(0.05);
      static_cast<D2DMatcher *>(matcher)->setNumLayers(1);
      static_cast<D2DMatcher *>(matcher)->setInputSource(source_grid);
      static_cast<D2DMatcher *>(matcher)->setInputTarget(target_grid);
    }
    if (mode_type == "corr") {
      matcher->setInputSource(source_grid->getMeans());
      matcher->setInputTarget(target_grid->getMeans());
    }
    if (mode_type == "robust") {
      static_cast<RobustMatcher *>(matcher)->setInputSource(source_grid);
      static_cast<RobustMatcher *>(matcher)->setInputTarget(target_grid);
    }

    if (mode_type == "corr_full") {
      matcher->setInputSource(scans[source_id]);
      matcher->setInputTarget(scans[target_id]);
    }
    if (mode_type == "basic") {
      // matcher->setInputSource(scans[source_id]);
      static_cast<P2DMatcher *>(matcher)->setNumLayers(1);
      matcher->setInputSource(source_grid->getMeans());
      static_cast<P2DMatcher *>(matcher)->setInputTarget(target_grid);
    }
    if (mode_type == "p2d") {
      matcher->setInputSource(scans[source_id]);
      static_cast<P2DMatcher *>(matcher)->setInputTarget(target_grid);
    }
    if (mode_type == "icpNDT") {
      matcher->setInputSource(source_grid->getMeans());
      matcher->setInputTarget(target_grid->getMeans());
    }
  }
  matcher->align(*output_m);
  e_match = std::chrono::system_clock::now();

  // std::cout << "sucess" << matcher->hasConverged() << std::endl;
  // std::cout << "result score: " << matcher->getFitnessScore() << std::endl;
  eigt::pose2d_t<double> calc_trans =
      eigt::getPoseFromTransform(eigt::convertToTransform<double>(
          matcher->getFinalTransformation().cast<double>()));

  elapsed_seconds = (e_match - s_match);
  // std::cout << "CALC TRANSFORM:" << calc_trans.transpose()
  //           << " calc time: " << elapsed_seconds.count() << std::endl;
  // std::cout << "DIFF:" << (real_trans - calc_trans).cwiseAbs().transpose()
  //           << std::endl
  //           << std::endl;
  total_error_ += (real_trans - calc_trans).cwiseAbs();
  // pcl::visualizePcl<pcl::PointXYZ>(scans[target_id], output_m);
}

void test(size_t start)
{
  // matches.clear();
  size_t target = 0;
  size_t iterations = 0;
  std::chrono::time_point<std::chrono::system_clock> s_total, e_total;
  std::chrono::duration<double> elapsed_seconds;
  size_t max_matchings = 1000;
  std::cout << "start test" << std::endl;
  s_total = std::chrono::system_clock::now();
  for (size_t i = start; i < scans.size(); ++i) {
    eigt::transform2d_t<double> real_trans =
        eigt::getTransFromPose(real_poses[target]).inverse() *
        eigt::getTransFromPose(real_poses[i]);

    if (std::abs(eigt::getAngle(real_trans)) > MIN_ROTATION &&
        eigt::getDisplacement(real_trans) > MIN_DISPLACEMENT) {
      // std::cout << "real trans: "
      //           << eigt::getPoseFromTransform(real_trans).transpose()
      //           << std::endl;
      // std::cout << std::abs(eigt::getAngle(real_trans)) << ":"
      //           << eigt::getDisplacement(real_trans) << std::endl;
      testMatch(i, target);
      target = i;
      ++iterations;
      if (iterations > max_matchings)
        break;
      // return;
    }
  }
  e_total = std::chrono::system_clock::now();
  elapsed_seconds = (e_total - s_total);
  Eigen::Vector3d error = total_error_ / iterations;

  std::cout << "-------------------------------------------------------------"
               "\n";
  std::cout << "average error: (x,y,theta)" << error.transpose()
            << " number of matchings:" << iterations
            << " average time per match: "
            << elapsed_seconds.count() / iterations << std::endl;
  std::cout << "-------------------------------------------------------------"
               "\n";
}

void testSkip(size_t start)
{
  size_t target = 0;
  std::cout << "start test" << std::endl;
  for (size_t i = start; i < scans.size(); i += 1) {
    eigt::transform2d_t<double> real_trans =
        eigt::transBtwPoses(real_poses[target], real_poses[i]);
    std::cout << "real trans: "
              << eigt::getPoseFromTransform(real_trans).transpose()
              << std::endl;
    std::cout << "matching " << i << "  " << target << std::endl;
    testMatch(i, target);
    target = i;
    // return;
  }
}

int main(int argc, char **argv)
{
  std::vector<std::string> args(argv, argv + argc);
  if (args.size() < 3 &&
      (args[2] != "d2d" && args[2] != "basic" && args[2] != "corr" &&
       args[2] != "corr_full" && args[2] != "robust" && args[2] != "proof" &&
       args[2] != "p2d" && args[2] != "icpNDT")) {
    std::cout << "Correct format of arguments: \n Path to the folder with "
                 "data.poses, data.scans was not provided\n Calculation engine "
                 "(d2d,basic,corr,robust) \n [Min displacement (m) min "
                 "rotation (rad)]";
    std::cout << std::endl;
    return 0;
  }
  if (args.size() == 4) {
    MIN_DISPLACEMENT = std::stod(args[3]);
  }
  if (args.size() == 5) {
    MIN_DISPLACEMENT = std::stod(args[3]);
    MIN_ROTATION = std::stod(args[4]);
  }
  // matcher->setResolution(1);
  preparePoseData(args[1]);
  prepareLaserScans(args[1]);
  // proofer =
  //     new D2DNormalDistributionsTransform2D<pcl::PointXYZ, pcl::PointXYZ>();
  proofer = new IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>();

  if (args[2] == "basic") {
    matcher = new P2DMatcher();
  } else if (args[2] == "d2d") {
    matcher = new D2DMatcher();
  } else if (args[2] == "corr") {
    matcher = new CorrMatcher();
    static_cast<CorrMatcher *>(matcher)->setCoarseStep(0.5);
  } else if (args[2] == "corr_full") {
    matcher = new CorrMatcher();
    static_cast<CorrMatcher *>(matcher)->setCoarseStep(0.5);
  } else if (args[2] == "robust") {
    matcher = new RobustMatcher();
  } else if (args[2] == "proof") {
    matcher = new IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>();
  } else if (args[2] == "p2d") {
    matcher = new P2DMatcher();
  } else if ((args[2] == "icpNDT")) {
    matcher = new IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>();
  }
  mode_type = args[2];
  // testMatch(4, 2);
  // testMatch(2, 0);
  // testMatch(6, 0);
  // testMatch(0, 10);
  // testMatch(0, 100);
  // size_t start = 80;
  test(0);
  delete matcher;
  delete proofer;
  return 0;
}
