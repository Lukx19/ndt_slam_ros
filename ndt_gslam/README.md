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
    
## NDT registration algorithms as a plugin to PointCloud library
### Basic usage
    #include <ndt_gslam/registration/d2d_ndt2d.h>

    int main(int argc, char ** argv)
    {
        typedef pcl::PointCloud<pcl::PointXYZ> Pcl;
        Pcl:Ptr target(new Pcl());
        Pcl:Ptr source(new Pcl());
        Pcl out_pcl;

        pcl::D2DNormalDistributionsTransform2D<pcl::PointXYZ, pcl::PointXYZ> matcher;
        matcher.setInputSource(source);
        matcher.setInputTarget(target);
        matcher.align(out_pcl);
        matcher.getFinalTransformation();
    }

### Usage with NDT frames
    #include <ndt_gslam/registration/d2d_ndt2d.h>
    #include <ndt_gslam/ndt/ndt_grid2d.h>
    #include <ndt_gslam/ndt/ndt_cell.h>

    int main(int argc, char ** argv)
    {
        typedef slamuk::NDTGrid2D<slamuk::NDTCell, pcl::PointXYZ> NDTFrame;
        NDTFrame::Ptr target(new NDTFrame());
        NDTFrame::Ptr source(new NDTFrame());
        pcl::PointCloud<pcl::PointXYZ> out_pcl;

        pcl::D2DNormalDistributionsTransform2D<pcl::PointXYZ, pcl::PointXYZ> matcher;
        matcher.setInputSource(source);
        matcher.setInputTarget(target);
        matcher.align(out_pcl);
        matcher.getFinalTransformation();
    }
The same API can be used for other registration algorithms in this package.
