# NDT-GSLAM
This is experimental version of NDT scanmatching frontend with Graph based backend SLAM.
## Requirements
- tested on ROS Kinetic
- g2o library:  apt-get install ros-kinetic-ligg2o
- Eigen3: usually included in ROS
- OpenMP compiler support
- OpenCV from official ubuntu repository is sufficient
- Pcl 1.7 (version from Ubuntu 14.04 (Trusty) repository is not working. PCL needs to be compiled with c++11 flags.)

## Installation
- Full desktop ROS installation usually needs these additional packages:
  - apt-get install ros-kinetic-tf ros-kinetic-laser-geometry ros-kinetic-pcl-ros tf-conversions
- catkin_make

## Quick Start
- download bag from MIT Stata dataset (http://projects.csail.mit.edu/stata/downloads.php) 
    - http://infinity.csail.mit.edu/data/2012/utilities/2012-01-28-11-12-01/2012-01-28-11-12-01.bag.noimages
- run bag in one terminal
- run ndt-glsm in second terminal
    - roslaunch ndt_gslam graph_slam_bag_MIT.launch
