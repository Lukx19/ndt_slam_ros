#include <laser_geometry/laser_geometry.h>
#include <ndt_gslam/ndt/cell_policy2d.h>
#include <ndt_gslam/ndt/ndt_cell.h>
#include <ndt_gslam/ndt/ndt_grid2d.h>
#include <ndt_gslam/ndt/ndt_mapper.h>
#include <ndt_gslam/ndt_slam_algorithm.h>
#include <ndt_gslam/registration/d2d_ndt2d.h>
#include <ndt_gslam/registration/ndt2d.h>
#include <pcl/common/time.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl_conversions/pcl_conversions.h>
#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <sensor_msgs/PointCloud.h>
#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <string>
#include <vector>

#include <ndt_gslam/utils/eigen_tools.h>
#include <ndt_gslam/utils/point_cloud_tools.h>
#include <ndt_gslam/utils/string_tools.h>

using namespace slamuk;

typedef pcl::PointCloud<pcl::PointXYZ> PointCloud;
typedef NDTGrid2D<NDTCell<CellPolicy2d>, pcl::PointXYZ> GridType;
typedef GridType::ConstPtr GridTypeConstPtr;
typedef GridType::Ptr GridTypePtr;
typedef pcl::D2DNormalDistributionsTransform2D<pcl::PointXYZ, pcl::PointXYZ>
    RobustMatcher;
typedef pcl::NormalDistributionsTransform2DEx<pcl::PointXYZ, pcl::PointXYZ>
    SimpleMatcher;
typedef NdtSlamAlgortihm SLAM;

double OUTLIER = 0.8;
size_t LAYERS = 4;
double STEP_SIZE = 0.01;
size_t SKIP = 3;

double EPSILON = 0.001;
double MIN_DISPLACEMENT = 0.1;
double MIN_ROTATION = 0.1;
size_t MAX_LINES = 5000;
float CELL_SIZE = 0.4;
bool initialized = false;
std::vector<typename PointCloud::Ptr> scans;
std::vector<Eigen::Vector3d> real_poses;
RobustMatcher matcher;
SimpleMatcher matcher_s;
GridTypePtr map;
eigt::transform2d_t<double> initial_trans;
Eigen::Matrix4f cummulative_trans;
bool flag = false;
SLAM algo;

Eigen::Vector3d old_origin;

GridTypePtr running_win;

size_t seq = 0;
PointCloud::Ptr all_points(new PointCloud());

void prepareLaserScans(std::string folder)
{
  sensor_msgs::PointCloud2 laser_pcl_msg;
  PointCloud laser_pcl;
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
  std::cout << "imported scans:" << scans.size() << std::endl;
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
  while (std::getline(pose_msg, line) && line_count < MAX_LINES * 2) {
    ++line_count;
    std::vector<std::string> parts = split(line, ",");
    pose << std::atof(parts[1].c_str()), std::atof(parts[2].c_str()),
        std::atof(parts[3].c_str());
    real_poses.push_back(pose);
  }
  pose_msg.close();
  std::cout << "imported poses:" << real_poses.size() << std::endl;
}
void visualizePcl(const PointCloud::ConstPtr &pcl)
{
  // Initializing point cloud visualizer
  boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(
      new pcl::visualization::PCLVisualizer("3D Viewer"));
  viewer->setBackgroundColor(0, 0, 0);

  pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> target_color(
      pcl, 255, 0, 0);
  viewer->addPointCloud<pcl::PointXYZ>(pcl, target_color, "vis cloud");
  viewer->setPointCloudRenderingProperties(
      pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "vis cloud");
  viewer->addCoordinateSystem(1.0, "global");
  // Wait until visualizer window is closed.
  while (!viewer->wasStopped()) {
    viewer->spinOnce(100);
    boost::this_thread::sleep(boost::posix_time::microseconds(100000));
  }
  viewer->close();
}

void visualizePcl(const PointCloud::ConstPtr &pcl1,
                  const PointCloud::ConstPtr &pcl2)
{
  // Initializing point cloud visualizer
  boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(
      new pcl::visualization::PCLVisualizer("3D Viewer"));
  viewer->setBackgroundColor(0, 0, 0);
  // visualize firs cloud
  pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> first_color(
      pcl1, 255, 0, 0);
  viewer->addPointCloud<pcl::PointXYZ>(pcl1, first_color, "first cloud");
  viewer->setPointCloudRenderingProperties(
      pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "first cloud");
  // visualize second cloud
  pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> second_color(
      pcl2, 0, 255, 0);
  viewer->addPointCloud<pcl::PointXYZ>(pcl2, second_color, "second cloud");
  viewer->setPointCloudRenderingProperties(
      pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "second cloud");
  viewer->addCoordinateSystem(1.0, "global");
  // Wait until visualizer window is closed.
  while (!viewer->wasStopped()) {
    viewer->spinOnce(100);
    boost::this_thread::sleep(boost::posix_time::microseconds(100000));
  }
  viewer->close();
}

void testMergeSimple(size_t scan_id)
{
  if (!initialized) {
    map.reset(new GridType());
    map->setCellSize(CELL_SIZE);
    map->setOrigin(Eigen::Vector3d(0, 0, 0));
    map->initialize(*scans[scan_id]);
    *all_points += *scans[scan_id];
    initial_trans =
        eigt::transBtwPoses(real_poses[scan_id], Eigen::Vector3d(0, 0, 0));
    cummulative_trans.setIdentity();
    initialized = true;
    std::cout << "INIT MAP: \n" << *map << std::endl << std::endl;
    return;
  }
  typename PointCloud::Ptr transformed_pcl(new PointCloud());
  GridTypePtr temp_grid(new GridType());
  temp_grid->setCellSize(CELL_SIZE);
  auto good_origin = eigt::transformPose(real_poses[scan_id], initial_trans);
  temp_grid->setOrigin(Eigen::Vector3d(0, 0, 0));
  temp_grid->mergeIn(*scans[scan_id], Eigen::Vector3d(0, 0, 0), true);
  matcher_s.setInputTarget(all_points);
  matcher_s.setInputSource(scans[scan_id]);
  // auto guess = eigt::transBtwFrames(good_origin, Eigen::Vector3d(0, 0, 0));
  matcher_s.align(*transformed_pcl, cummulative_trans);
  cummulative_trans = matcher.getFinalTransformation();
  pcl::transformPointCloud(*scans[scan_id], *transformed_pcl,
                           cummulative_trans);

  map->mergeInTraced(*transformed_pcl, Eigen::Vector3d(0, 0, 0), true);

  // std::cout << "MAP AFTER ADDITION: " << seq << " \n"
  //           << *map << std::endl
  //           << std::endl;
  // visualizePcl(all_points, transformed_pcl);
  *all_points += *transformed_pcl;
  seq++;
}

void testMergeD2D(size_t scan_id)
{
  if (!initialized) {
    map.reset(new GridType());
    map->setCellSize(CELL_SIZE);
    map->setOrigin(Eigen::Vector3d(0, 0, 0));
    map->initialize(*scans[scan_id]);
    *all_points += *scans[scan_id];
    initial_trans =
        eigt::transBtwPoses(real_poses[scan_id], Eigen::Vector3d(0, 0, 0));
    cummulative_trans.setIdentity();
    initialized = true;
    std::cout << "INIT MAP: \n" << *map << std::endl << std::endl;
    return;
  }
  typename PointCloud::Ptr transformed_pcl(new PointCloud());
  GridTypePtr temp_grid(new GridType());
  temp_grid->setCellSize(CELL_SIZE);
  auto good_origin = eigt::transformPose(real_poses[scan_id], initial_trans);
  temp_grid->setOrigin(Eigen::Vector3d(0, 0, 0));
  temp_grid->mergeIn(*scans[scan_id], Eigen::Vector3d(0, 0, 0), true);
  matcher.setInputTarget(map);
  matcher.setInputSource(temp_grid);
  // auto guess = eigt::transBtwPoses(good_origin, Eigen::Vector3d(0, 0, 0));
  matcher.align(*transformed_pcl, cummulative_trans);
  cummulative_trans = matcher.getFinalTransformation();
  pcl::transformPointCloud(*scans[scan_id], *transformed_pcl,
                           cummulative_trans);

  map->mergeInTraced(*transformed_pcl, Eigen::Vector3d(0, 0, 0), true);

  // std::cout << "MAP AFTER ADDITION: " << seq << " \n"
  //           << *map << std::endl
  //           << std::endl;
  // visualizePcl(all_points, transformed_pcl);
  *all_points += *transformed_pcl;
  seq++;
}

void runningWindow(size_t scan_id)
{
  if (!initialized) {
    map.reset(new GridType());
    map->setCellSize(CELL_SIZE);
    map->setOrigin(Eigen::Vector3d(0, 0, 0));
    map->initialize(*scans[scan_id]);
    running_win->setCellSize(CELL_SIZE);
    running_win->initialize(*scans[scan_id]);
    *all_points += *scans[scan_id];
    initial_trans =
        eigt::transBtwPoses(real_poses[scan_id], Eigen::Vector3d(0, 0, 0));
    cummulative_trans.setIdentity();
    initialized = true;
    std::cout << "INIT MAP: \n" << *map << std::endl << std::endl;
    seq = 0;
    old_origin.setZero();
    matcher.setOulierRatio(OUTLIER);
    matcher.setStepSize(STEP_SIZE);
    matcher.setCellSize(0.4);
    matcher.setNumLayers(LAYERS);

    matcher_s.setOulierRatio(OUTLIER);
    matcher_s.setStepSize(STEP_SIZE);
    matcher_s.setNumLayers(LAYERS);
    matcher_s.setCellSize(0.4);
    return;
  }
  typename PointCloud::Ptr transformed_pcl(new PointCloud());

  GridTypePtr temp_grid(new GridType());
  temp_grid->setCellSize(CELL_SIZE);

  temp_grid->setOrigin(running_win->getOrigin());

  temp_grid->mergeIn(*scans[scan_id], running_win->getOrigin(), true);
  // matcher.setInputTarget(running_win);
  // matcher.setInputSource(temp_grid);
  // matcher.align(*transformed_pcl, cummulative_trans);

  // cummulative_trans = matcher.getFinalTransformation();

  matcher_s.setInputTarget(running_win);
  matcher_s.setInputSource(temp_grid->getMeans());
  matcher_s.align(*transformed_pcl, cummulative_trans);

  cummulative_trans = matcher_s.getFinalTransformation();

  // if (flag == true) {
  //   flag = false;
  //   visualizePcl(running_win->getMeans(), transformed_pcl);
  // }

  eigt::transform2d_t<double> cumul_trans2d =
      eigt::convertToTransform<float>(cummulative_trans).cast<double>();

  pcl::transformPointCloud(*scans[scan_id], *transformed_pcl,
                           eigt::convertFromTransform(cumul_trans2d));
  temp_grid->transform(cumul_trans2d);
  running_win->mergeInTraced(*temp_grid, true, false);
  // running_win->mergeInTraced(*transformed_pcl, running_win->getOrigin());
  map->mergeInTraced(*temp_grid, true, true);

  cumul_trans2d = running_win->move(cumul_trans2d);
  cummulative_trans = eigt::convertFromTransform(cumul_trans2d).cast<float>();

  // std::cout << "running window seq:" << seq << std::endl
  //           << *running_win << std::endl;

  //*all_points += *temp_grid->getMeansTransformed();
  // if (std::abs(old_origin.sum() - running_win->getOrigin().sum()) > 0.1) {
  //   visualizePcl(map->getMeans(), running_win->getMeansTransformed());
  //   old_origin = running_win->getOrigin();
  //   flag = true;
  // }
  seq++;
}

void moveWindow()
{
  GridTypePtr temp_grid(new GridType());
  temp_grid->setCellSize(CELL_SIZE);

  temp_grid->setOrigin(Eigen::Vector3d(0, 0, 0));

  temp_grid->mergeIn(*scans[80], Eigen::Vector3d(0, 0, 0), true);
  eigt::transform2d_t<double> trans;
  trans.matrix().block(0, 2, 2, 1) << 4, 0;
  std::cout << trans.matrix();
  std::cout << *temp_grid << std::endl;
  temp_grid->move(trans);
  std::cout << *temp_grid << std::endl;
  visualizePcl(scans[80], temp_grid->getMeansTransformed());
  visualizePcl(scans[80], temp_grid->getMeans());
}

void transformGrid()
{
  GridTypePtr temp_grid(new GridType());
  temp_grid->setCellSize(CELL_SIZE);

  temp_grid->setOrigin(Eigen::Vector3d(1, 3, 3.12));

  temp_grid->mergeIn(*scans[80], Eigen::Vector3d(0, 0, 0), true);
  eigt::transform2d_t<double> trans =
      eigt::getTransFromPose(Eigen::Vector3d(3, 0, 1.58));

  std::cout << trans.matrix();
  std::cout << *temp_grid << std::endl;
  temp_grid->transform(trans);
  std::cout << *temp_grid << std::endl;
  visualizePcl(scans[80], temp_grid->getMeansTransformed());
  visualizePcl(scans[80], temp_grid->getMeans());
}

void runSLAM()
{
  size_t target = 0;
  algo.update(SLAM::Transform::Identity(), SLAM::Covar::Identity(), *scans[0],
              ros::Time(static_cast<double>(0)));
  for (size_t i = 1; i < scans.size() - 1; ++i) {
    eigt::transform2d_t<double> odom_trans =
        eigt::transBtwPoses(real_poses[target], real_poses[i]);
    algo.update(odom_trans, SLAM::Covar::Identity(), *scans[i],
                ros::Time(static_cast<double>(i)));
    algo.getOccupancyGrid("aaaa");
    algo.getPclMap("aaa");
    algo.getPclMap2("aaa");
    target = i;
  }
}

void mappingKnownPoses()
{
  GridTypePtr merge_map(new GridType());
  merge_map->setCellSize(0.25);
  std::cout << "start test" << std::endl;
  for (size_t i = 7000; i < scans.size() - 10000; i += 1) {
    Eigen::Vector3d glob_pose = eigt::getPoseFromTransform(
        eigt::getTransFromPose(real_poses[0]).inverse() *
        eigt::getTransFromPose(real_poses[i]));
    merge_map->mergeInTraced(*scans[i], glob_pose, true);
  }
  std::cout << *merge_map << std::endl;
  pcl::visualizePcl<pcl::PointXYZ>(merge_map->getMeansTransformed());
}

void test(size_t start)
{
  {
    pcl::ScopeTime t_init("[NDT_GRID2D_MERGING]: total calculation time:");
    // matches.clear();
    size_t target = start;
    runningWindow(target);
    // testMergeD2D(start);
    std::cout << "start test" << std::endl;
    {
      pcl::ScopeTime t_recalc("calc time");
      for (size_t i = start; i < scans.size(); i += SKIP) {
        eigt::transform2d_t<double> real_trans =
            eigt::transBtwPoses(real_poses[target], real_poses[i]);
        // if (eigt::getAngle(real_trans) > MIN_ROTATION ||
        //     eigt::getDisplacement(real_trans) > MIN_DISPLACEMENT) {
        // std::cout << "displace: " << eigt::getDisplacement(real_trans)
        //           << "rot: " << eigt::getAngle(real_trans) << std::endl;
        // testMergeD2D(i);
        target = i;
        runningWindow(i);
        // }
      }
    }
  }
  std::cout << "Running win: " << seq << " \n"
            << *running_win << std::endl
            << std::endl;
  std::cout << "MAP AFTER ADDITION: " << seq << " \n"
            << *map << std::endl
            << std::endl;

  visualizePcl(map->getMeansTransformed(), running_win->getMeansTransformed());
}

void testGridMerginTransformation()
{
  size_t scan_id = 80;
  GridTypePtr temp_grid(new GridType());
  temp_grid->setCellSize(CELL_SIZE);
  temp_grid->setOrigin(Eigen::Vector3d(10, 10, 3.12));
  // scan is captured in robots frame which is always in 0,0,0
  temp_grid->mergeInTraced(*scans[scan_id], Eigen::Vector3d(0, 0, 0), true);
  std::cout << *temp_grid << std::endl;
  visualizePcl(scans[scan_id], temp_grid->getMeansTransformed());
}

void testGridMerginTransformation2()
{
  size_t scan_id = 80;
  GridTypePtr temp_grid(new GridType());
  temp_grid->setCellSize(CELL_SIZE);
  temp_grid->setOrigin(Eigen::Vector3d(10, 10, 3.12));
  // scan is captured in robots frame which is always in 0,0,0
  temp_grid->mergeInTraced(*scans[scan_id], Eigen::Vector3d(10, 10, 3.12),
                           true);
  std::cout << *temp_grid << std::endl;
  visualizePcl(scans[scan_id], temp_grid->getMeansTransformed());
}

void testGridMerginFrames(bool traced, bool transform)
{
  size_t scan_id = 80;
  GridTypePtr temp_grid(new GridType());
  GridTypePtr base_grid(new GridType());
  temp_grid->setCellSize(CELL_SIZE);
  base_grid->setCellSize(CELL_SIZE);

  temp_grid->setOrigin(Eigen::Vector3d(10, 10, 3.12));
  base_grid->setOrigin(Eigen::Vector3d(-10, -10, 0));

  temp_grid->mergeIn(*scans[scan_id], Eigen::Vector3d(10, 10, 3.12), true);
  base_grid->mergeIn(*scans[scan_id], Eigen::Vector3d(-10, -10, 0), true);
  std::cout << *temp_grid << std::endl << std::endl;
  std::cout << *base_grid << std::endl << std::endl;
  if (traced)
    base_grid->mergeInTraced(*temp_grid, transform, true);
  else
    base_grid->mergeIn(*temp_grid, transform, true);

  base_grid->mergeIn(*temp_grid, false, true);
  std::cout << *base_grid << std::endl << std::endl;
  visualizePcl(base_grid->getMeansTransformed(),
               temp_grid->getMeansTransformed());
}

void testDynamicEnviro()
{
  size_t scan_id = 80;
  GridTypePtr temp_grid(new GridType());
  GridTypePtr base_grid(new GridType());
  GridTypePtr valid_grid(new GridType());

  temp_grid->setCellSize(CELL_SIZE);
  base_grid->setCellSize(CELL_SIZE);
  valid_grid->setCellSize(CELL_SIZE);

  temp_grid->setOrigin(Eigen::Vector3d(0, 0, 0));
  base_grid->setOrigin(Eigen::Vector3d(0, 0, 0));
  valid_grid->setOrigin(Eigen::Vector3d(0, 0, 0));

  temp_grid->mergeIn(*scans[scan_id], Eigen::Vector3d(-1, 0, 0), true);
  base_grid->mergeIn(*scans[scan_id], Eigen::Vector3d(0, 0, 0), true);
  valid_grid->mergeIn(*scans[scan_id], Eigen::Vector3d(0, 0, 0), true);

  std::cout << *temp_grid << std::endl << std::endl;
  std::cout << *base_grid << std::endl << std::endl;
  base_grid->mergeInTraced(*temp_grid, true);
  std::cout << *base_grid << std::endl << std::endl;
  for (size_t i = 0; i < 1; ++i) {
    base_grid->mergeInTraced(*valid_grid, true);
    std::cout << *base_grid << std::endl << std::endl;
  }
}

void testMapMerging()
{
  NDTMapper<NDTCell<CellPolicy2d>, pcl::PointXYZ> merged_map;
  std::vector<GridTypePtr> grids;
  grids.emplace_back(new GridType(Eigen::Vector3d(0, 0, 2)));
  grids.emplace_back(new GridType(Eigen::Vector3d(0, 4, 0)));
  grids.emplace_back(new GridType(Eigen::Vector3d(0, -4, 0)));
  grids.emplace_back(new GridType(Eigen::Vector3d(0, 0, 3.14)));
  grids.emplace_back(new GridType(Eigen::Vector3d(-6, 0, 3.14)));
  double t = 0;
  for (auto &frame : grids) {
    frame->mergeIn(*scans[80], frame->getOrigin());
    merged_map.addFrame(frame, ros::Time(t));
    ++t;
  }
  merged_map.recalc(ros::Time(100));
  // cv::namedWindow("MAP", cv::WINDOW_NORMAL);
  // cv::imshow("MAP", merged_map.getOccupancyMap());
  // cv::waitKey(0);
  for (size_t i = 0; i < 300; ++i) {
    GridTypePtr temp_grid(new GridType(Eigen::Vector3d(0, 0, 0)));
    temp_grid->mergeInTraced(*scans[80], temp_grid->getOrigin());
    merged_map.addFrame(temp_grid, ros::Time((double)i));
  }
  {
    pcl::ScopeTime t_recalc("recalc:time 100");
    merged_map.recalc(ros::Time(1000));
  }
}

int main(int argc, char **argv)
{
  std::vector<std::string> args(argv, argv + argc);
  if (args.size() < 2) {
    std::cout << "Correct format of arguments: \n Path to the folder with "
                 "data.poses, data.scans was not provided";
    std::cout << std::endl;
    return 0;
  }
  if (args.size() >= 3)
    MAX_LINES = std::stoi(args[2]);
  float win_size = 15;
  if (args.size() >= 4)
    win_size = std::stof(args[3]);
  if (args.size() >= 8) {
    OUTLIER = std::stod(args[4]);
    STEP_SIZE = std::stod(args[5]);
    SKIP = std::stoi(args[6]);
    LAYERS = std::stoi(args[7]);
  }
  preparePoseData(args[1]);
  prepareLaserScans(args[1]);
  running_win = GridTypePtr(new GridType());
  running_win->setCellSize(CELL_SIZE);
  running_win->setOrigin(Eigen::Vector3d(0, 0, 0));
  running_win->enlarge(-win_size, -win_size, win_size, win_size);
  // testGridMerginTransformation();
  // testGridMerginTransformation2();
  // testGridMerginFrames(false, false);
  // testGridMerginFrames(false, true);
  // testGridMerginFrames(true, false);
  // testGridMerginFrames(true, true);
  // testDynamicEnviro();
  // testMapMerging();
  // moveWindow();
  // transformGrid();
  // size_t start = 200;
  // test(0);
  // runSLAM();
  // testLoopTools();
  mappingKnownPoses();
  return 0;
}